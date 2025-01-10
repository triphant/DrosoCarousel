classdef CaroDataLoader < BasicDataHandleList & CaroInterface
%CAROINTERFACE Class to load carousel data, repair problems and return the
%   processed data
%
%   Author: Tilman Triphan, tilman.triphan@uni-leipzig.de
%   License: GNU General Public License v3.0
%   Copyright (C) 2017-2025  Tilman Triphan
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.

    properties
        proc
    end

    properties (Constant = true, Hidden = true)
        cdl_version_code = 'cdl_20230606';
        CAROTYPES = {'singleCaro','doubleCaro'}
        DATATYPES = {'.csv','.h5'}
        MAXFLYLENGTH = 20;
    end

    methods

        function obj = CaroDataLoader(data)
            if nargin > 0
                obj.data = data;
            end
            obj.cdl_version_data = obj.cdl_version_code;
            obj.proc.dataLoaded = false;
            obj.proc.loading = {};
            obj.proc.errorhandling = {};
        end

        function loadMetaData(obj,filename)
            % loadMetaData(obj,filename) loads a metadata file (.csv or 
            % .json)
            arguments
                obj (1,1)
                filename char {mustBeFile} = obj.info.metaDataFile;
            end
            if endsWith(filename,'csv')
                mdt = importdata(filename);
                mds = cell2struct(strsplit(cell2mat(mdt(2,:)),','),strsplit(cell2mat(mdt(1,:)),','),2);
                obj.info.type = 'SingleCaro';
                obj.info.expTime = datetime([mds.date,mds.time],'InputFormat','yy/MM/ddHH:mm');
                obj.info.rigID = mds.Rig;
                obj.info.fps = single(str2double(mds.FPS));
                obj.info.genotype = mds.Genotype_m;
                obj.info.age = single(str2double(mds.Age_m_days) + str2double(mds.Age_m_hours)/24);
                obj.info.temperature = single(str2double(mds.Temp_degC));
                obj.info.humidity = single(str2double(mds.Humidity_perc));
                obj.info.motorOnOff = logical(str2double(mds.Motor_OnOff));
                obj.info.motorDir = char(string(mds.Motor_direction).replace('RIGHT','cw').replace('LEFT','ccw').replace('cc','ccw'));
                tempStr = string(mds.Add_Information);
                obj.info.diskPos = tempStr.contains('higher')-tempStr.contains('lower')-tempStr.contains('tiefer');
                obj.info.food = mds.Food;

                obj.ctrl.maxFrame = single(str2double(mds.Videolength_sec)*obj.info.fps);
            elseif endsWith(filename,'json')
                j = jsondecode(fileread(filename));
                obj.info.type = 'DoubleCaro';
                obj.info.expTime = datetime([j.date,j.time],'InputFormat','yyMMddHHmmss');
                obj.info.rigID = j.rig;
                obj.info.fps = single(j.FPS);
                obj.info.genotype = j.genotype;
                obj.info.age = single(j.age);
                obj.info.temperature = single(j.temperature);
                obj.info.humidity = single(j.humidity);
                obj.info.motorDir1 = j.motor1_direction;
                obj.info.motorDir2 = j.motor2_direction;
                obj.info.food = j.food;

                obj.ctrl.maxFrame = single(j.video_length);
            else
                disp('invalid metadata file!')
                return
            end

            obj.ctrl.framesPerSec = obj.info.fps;
            obj.ctrl.videoLength = obj.ctrl.maxFrame / 60 / obj.info.fps;
            obj.ctrl.imageResX = 768;
            obj.ctrl.imageResY = 768;
            obj.ctrl.mm2pix = 5.9881;

            fn = cell2mat(obj.getInfo('fullName'));
            switch obj.info.type
                case 'SingleCaro'
                    obj.info.name = fn(1:14);
                case 'DoubleCaro'
                    obj.info.name = fn(1:17);
            end
        end

        function setTimeAndDate(obj)
            % setTimeAndDate(obj) sets date and time
            obj.info.startTime = obj.info.expTime;
            vl = obj.ctrl.videoLength;
            if obj.info.startTime.Minute < vl
                obj.info.startTime.Minute = vl;
            else
                obj.info.startTime.Hour = obj.info.startTime.Hour + 1;
                obj.info.startTime.Minute = 0;
            end
            obj.info.startTime.Second = 0;
            st = obj.info.startTime;
            obj.addFieldAndValues('dateTime', ...
                st:minutes(vl):st+minutes(vl)*(obj.getLength-1));
        end

        function setRoi(obj,roiname,x,y,rmin,rmax,roitype)
            % setRoi(obj,roiname,x,y,rmin,rmax,roitype) updates ROIs
            obj.ctrl.rois.(roiname).x = x;
            obj.ctrl.rois.(roiname).y = y;
            obj.ctrl.rois.(roiname).rmin = rmin;
            obj.ctrl.rois.(roiname).rmax = rmax;
            obj.ctrl.rois.(roiname).type = roitype;
        end

        function loadRoiSet(obj,filename)
            % loadRoiSet(obj,filename) loads roi definitions and creates
            % ROIs
            arguments
                obj (1,1)
                filename char {mustBeFile}
            end
            rois = readtable(filename,'FileType','delimitedtext');
            switch cell2mat(obj.getInfo('type'))
                case 'SingleCaro'
                    x = ceil(rois.X(1));
                    y = ceil(rois.Y(1));
                    r = ceil(rois.Width(1)/2);
                    obj.setRoi('arena',x,y,0,r,'circ');
                    obj.setRoi('border',x,y,r-obj.conf.borderWidth,r-obj.conf.borderBuffer,'circ');
                    x = ceil(rois.X(2));
                    y = ceil(rois.Y(2));
                    r = ceil(rois.Width(2)/2);
                    obj.setRoi('disk',x,y,0,r,'circ');
                    obj.setRoi('ring',x,y,r,r+obj.conf.ringWidth,'circ');
                    x = ceil(rois.X(3));
                    y = ceil(rois.Y(3));
                    r = ceil(rois.Width(3)/2);
                    obj.setRoi('food',x,y,0,r,'circ');
                case 'DoubleCaro'
                    x = ceil(rois.X(1));
                    y = ceil(rois.Y(1));
                    r = ceil(rois.Width(1)/2);
                    obj.setRoi('arena',x,y,0,r,'circ');
                    obj.setRoi('border',x,y,r-obj.conf.borderWidth,r-obj.conf.borderBuffer,'circ');
                    x = ceil(rois.X(2));
                    y = ceil(rois.Y(2));
                    r = ceil(rois.Width(2)/2);
                    obj.setRoi('disk1',x,y,0,r,'circ');
                    obj.setRoi('ring1',x,y,r,r+obj.conf.ringWidth,'circ');
                    x = ceil(rois.X(3));
                    y = ceil(rois.Y(3));
                    r = ceil(rois.Width(3)/2);
                    obj.setRoi('disk2',x,y,0,r,'circ');
                    obj.setRoi('ring2',x,y,r,r+obj.conf.ringWidth,'circ');
                    x = ceil(rois.X(4));
                    y = ceil(rois.Y(4));
                    r = ceil(rois.Width(4)/2);
                    obj.setRoi('food',x,y,0,r,'circ');
            end
        end

        function setRoiDefaultValues(obj,expType)
            % setRoiDefaultValues(obj,expType) sets default values for ROIs
            obj.ctrl.rois.arena.color = [0.8 0.8 0.8];
            obj.ctrl.rois.arena.comment = 'whole arena';
            obj.ctrl.rois.arena.prio = -100;
            obj.ctrl.rois.arena.sort = 0;
            obj.ctrl.rois.arena.type = 'circ';
            obj.ctrl.rois.arena.use = false;

            obj.ctrl.rois.empty.color = [0.8 0.8 0.8];
            obj.ctrl.rois.empty.comment = 'arena excluding all defined rois';
            obj.ctrl.rois.empty.inc = {};
            obj.ctrl.rois.empty.prio = -100;
            obj.ctrl.rois.empty.rule = 'not';
            obj.ctrl.rois.empty.sort = 0;
            obj.ctrl.rois.empty.type = 'comb';
            obj.ctrl.rois.empty.use = false;

            obj.ctrl.rois.nan.color = [0 0 0];
            obj.ctrl.rois.nan.comment = 'invalid values';
            obj.ctrl.rois.nan.func = 'isNaN';
            obj.ctrl.rois.nan.prio = 0;
            obj.ctrl.rois.nan.sort = 5;
            obj.ctrl.rois.nan.type = 'func';
            obj.ctrl.rois.nan.use = true;

            obj.ctrl.rois.food.color = [0 0.7 0];
            obj.ctrl.rois.food.comment = 'food patch';
            obj.ctrl.rois.food.prio = 0;
            obj.ctrl.rois.food.sort = 2;
            obj.ctrl.rois.food.type = 'circ';
            obj.ctrl.rois.food.use = true;

            obj.ctrl.rois.border.color = [0 0.5 0.5];
            obj.ctrl.rois.border.comment = 'border around arena';
            obj.ctrl.rois.border.prio = -10;
            obj.ctrl.rois.border.sort = 4;
            obj.ctrl.rois.border.type = 'circ';
            obj.ctrl.rois.border.use = true;

            obj.ctrl.rois.nodata.color = [1 1 1];
            obj.ctrl.rois.nodata.comment = 'no data e.g. outside of experiment time window';
            obj.ctrl.rois.nodata.type = 'spec';
            obj.ctrl.rois.nodata.use = false;

            switch expType
                case 'SingleCaro'
                    obj.ctrl.rois.empty.inc = {
                        'arena','disk','food','water','ring','border'};

                    obj.ctrl.rois.disk.color = [1 0 0];
                    obj.ctrl.rois.disk.comment = 'carousel disk';
                    obj.ctrl.rois.disk.prio = 10;
                    obj.ctrl.rois.disk.sort = 1;
                    obj.ctrl.rois.disk.type = 'circ';
                    obj.ctrl.rois.disk.use = true;

                    obj.ctrl.rois.ring.color = [0.5 0.5 0];
                    obj.ctrl.rois.ring.comment = 'ring around carousel disk';
                    obj.ctrl.rois.ring.prio = 0;
                    obj.ctrl.rois.ring.sort = 3;
                    obj.ctrl.rois.ring.type = 'circ';
                    obj.ctrl.rois.ring.use = true;

                case 'DoubleCaro'
                    obj.ctrl.rois.empty.inc = {
                        'arena','disk1','disk2','food','water','ring1','ring2','border'};

                    obj.ctrl.rois.disk1.color = [0.7 0 0];
                    obj.ctrl.rois.disk1.comment = 'left carousel disk';
                    obj.ctrl.rois.disk1.prio = 10;
                    obj.ctrl.rois.disk1.sort = 1;
                    obj.ctrl.rois.disk1.type = 'circ';
                    obj.ctrl.rois.disk1.use = true;

                    obj.ctrl.rois.disk2.color = [1 0 0];
                    obj.ctrl.rois.disk2.comment = 'right carousel disk';
                    obj.ctrl.rois.disk2.prio = 10;
                    obj.ctrl.rois.disk2.sort = 1;
                    obj.ctrl.rois.disk2.type = 'circ';
                    obj.ctrl.rois.disk2.use = true;

                    obj.ctrl.rois.ring1.color = [0.5 0.5 0];
                    obj.ctrl.rois.ring1.comment = 'ring around left carousel disk';
                    obj.ctrl.rois.ring1.prio = 0;
                    obj.ctrl.rois.ring1.sort = 3;
                    obj.ctrl.rois.ring1.type = 'circ';
                    obj.ctrl.rois.ring1.use = true;

                    obj.ctrl.rois.ring2.color = [0.8 0.8 0];
                    obj.ctrl.rois.ring2.comment = 'ring around right carousel disk';
                    obj.ctrl.rois.ring2.prio = 0;
                    obj.ctrl.rois.ring2.sort = 3;
                    obj.ctrl.rois.ring2.type = 'circ';
                    obj.ctrl.rois.ring2.use = true;

                    obj.ctrl.rois.turn1.color = [0 0 0];
                    obj.ctrl.rois.turn1.comment = 'left carousel is turning';
                    obj.ctrl.rois.turn1.func = 'isDiskMoving1';
                    obj.ctrl.rois.turn1.prio = 0;
                    obj.ctrl.rois.turn1.sort = 5;
                    obj.ctrl.rois.turn1.type = 'func';
                    obj.ctrl.rois.turn1.use = false;

                    obj.ctrl.rois.turn2.color = [0 0 0];
                    obj.ctrl.rois.turn2.comment = 'right carousel is turning';
                    obj.ctrl.rois.turn2.func = 'isDiskMoving2';
                    obj.ctrl.rois.turn2.prio = 0;
                    obj.ctrl.rois.turn2.sort = 5;
                    obj.ctrl.rois.turn2.type = 'func';
                    obj.ctrl.rois.turn2.use = false;

                    obj.ctrl.rois.diskturn1.color = [0.7 0 0];
                    obj.ctrl.rois.diskturn1.comment = 'fly on left carousel (turning)';
                    obj.ctrl.rois.diskturn1.inc = {'disk1','turn1'};
                    obj.ctrl.rois.diskturn1.prio = 50;
                    obj.ctrl.rois.diskturn1.rule = 'and';
                    obj.ctrl.rois.diskturn1.sort = 5;
                    obj.ctrl.rois.diskturn1.type = 'comb';
                    obj.ctrl.rois.diskturn1.use = false;

                    obj.ctrl.rois.diskturn2.color = [1 0 0];
                    obj.ctrl.rois.diskturn2.comment = 'fly on right carousel (turning)';
                    obj.ctrl.rois.diskturn2.inc = {'disk2','turn2'};
                    obj.ctrl.rois.diskturn2.prio = 50;
                    obj.ctrl.rois.diskturn2.rule = 'and';
                    obj.ctrl.rois.diskturn2.sort = 5;
                    obj.ctrl.rois.diskturn2.type = 'comb';
                    obj.ctrl.rois.diskturn2.use = false;

                    obj.ctrl.rois.diskstop1.color = [0 0 0.7];
                    obj.ctrl.rois.diskstop1.comment = 'fly on left carousel (stopped)';
                    obj.ctrl.rois.diskstop1.inc = {'disk1','turn1'};
                    obj.ctrl.rois.diskstop1.prio = -50;
                    obj.ctrl.rois.diskstop1.rule = 'not';
                    obj.ctrl.rois.diskstop1.sort = 5;
                    obj.ctrl.rois.diskstop1.type = 'comb';
                    obj.ctrl.rois.diskstop1.use = false;

                    obj.ctrl.rois.diskstop2.color = [0 0 1];
                    obj.ctrl.rois.diskstop2.comment = 'fly on right carousel (stopped)';
                    obj.ctrl.rois.diskstop2.inc = {'disk2','turn2'};
                    obj.ctrl.rois.diskstop2.prio = -50;
                    obj.ctrl.rois.diskstop2.rule = 'not';
                    obj.ctrl.rois.diskstop2.sort = 5;
                    obj.ctrl.rois.diskstop2.type = 'comb';
                    obj.ctrl.rois.diskstop2.use = false;
            end
        end

        function [cdl, errorlog] = loadData(obj,params)
            % [cdl, errorlog] = loadData(obj,params) loads the data
            arguments
                obj (1,1)
                params.maxFiles (1,1) = NaN
                params.minFiles (1,1) = 48
                params.preferType (1,:) char {mustBeMember(params.preferType,{'.h5','.csv'})} = '.h5'
            end

            obj.checkForError('mixedExperimentTypes',{1})
            errorlog = {};

            invalidData = false(numel(obj),1);
            nFiles = cell2mat(obj.getInfo('h5Count'));

            for o=1:numel(obj)
                if nFiles(o) < params.minFiles
                    invalidData(o) = 1;
                    errorlog(end+1) = {['Invalid length in ' cell2mat(obj(o).getInfo('fullName'))]}; %#ok<AGROW>
                    disp('Invalid length, skipping...')
                    continue
                end
                if strcmp(params.preferType,'.h5')
                    if nFiles(o) == 0
                        invalidData(o) = 1;
                        params.preferType = '.csv';
                        errorlog(end+1) = {['No h5 files in ' cell2mat(obj(o).getInfo('fullName'))]}; %#ok<AGROW>
                        disp('No h5 files...')
                        continue
                    end
                end
                datapath = cell2mat(obj(o).getInfo('dataPath'));
                files = dir(strcat(datapath,'\*',params.preferType));
                fprintf('Found %d files in folder %s, loading now...\n',numel(files),datapath)
                data = struct([]);
                fileSize = zeros(numel(files),1);

                for i=1:min(numel(files),params.maxFiles)
                    fprintf('Loading file %s (%d bytes)\n',files(i).name,files(i).bytes)
                    fileSize(i) = files(i).bytes;
                    h5struct = h5read(strcat(datapath,'\',files(i).name),'/df_with_missing/table');
                    h5data = h5struct.values_block_0';
                    % check if the experiment was tracked with head and
                    % tail or center of gravity settings
                    switch size(h5data,2)
                        case 3
                            data(i).x = single(h5data(:,1));
                            data(i).y = single(h5data(:,2));
                            data(i).qual = single(h5data(:,3));
                        case 6
                            data(i).xvals = single(h5data(:,[1 4]));
                            data(i).yvals = single(h5data(:,[2 5]));
                            data(i).qvals = single(h5data(:,[3 6]));
                            data(i).x = mean(data(i).xvals,2);
                            data(i).y = mean(data(i).yvals,2);
                            data(i).qual = mean(data(i).qvals,2);
                    end
                    data(i).frameCount = uint32(max(h5struct.index+1));
                    data(i).file = files(i).name;
                end
                obj(o).data = data;
                if isfield(obj(o).data(1),'xvals')
                    obj(o).ctrl.tracking = 'hat';
                    disp('dataset uses head and tail (hat) tracking')
                else
                    obj(o).ctrl.tracking = 'cog';
                    disp('dataset uses center of gravity (cog) tracking')
                end
                expType = cell2mat(obj(o).getInfo('type'));
                switch expType
                    case 'SingleCaro'
                        obj(o).loadRoiSet([datapath '\Roiset.xls']);
                        ids = cellfun(@(x) x(1:regexp(x,'\_[LD]_')+5), obj(o).get('file'),'UniformOutput', false);
                    case 'DoubleCaro'
                        obj(o).loadRoiSet([datapath '\RoisetDoppelcaro.xls']);
                        ids = cellfun(@(x) x(1:regexp(x,'DLC')-1), obj(o).get('file'),'UniformOutput', false);
                end
                obj(o).addFieldAndValues('fileID',ids);
                obj(o).setRoiDefaultValues(expType);
                obj(o).updateRoiConf();
                obj(o).setTimeAndDate();
                obj(o).info.tracking = obj(o).ctrl.tracking;
                obj(o).proc.dataLoaded = true;
                obj(o).proc.loading{end+1} = 'loadData';
            end
            cdl = obj(~invalidData);
        end

        function cdhl = toCaroDataHandleList(obj)
            % cdhl = toCaroDataHandleList(obj) converts the CaroDataLoader
            % to a CaroDataHandleList or DoubleCaroDataHandleList
            arguments
                obj (1,:)
            end
            obj.checkForError('mixedExperimentTypes',{1})
            switch cell2mat(obj(1).getInfo('type'))
                case 'SingleCaro'
                    className = 'CaroDataHandleList';
                case 'DoubleCaro'
                    className = 'DoubleCaroDataHandleList';
            end
            cdhl = feval([className,'.empty'],numel(obj),0);
            invalidData = false(numel(obj),1);
            for o=1:numel(obj)
                if isempty(obj(o).getData)
                    invalidData(o) = 1;
                    continue
                end
                cdhl(o) = feval(className,obj(o).getData); %#ok<FVAL>
                cdhl(o).info = obj(o).info;
                cdhl(o).conf = obj(o).conf;
                cdhl(o).ctrl = obj(o).ctrl;
                cdhl(o).proc = obj(o).proc;
                cdhl(o).conf.indexMode = 'video';
            end
            cdhl = cdhl(~invalidData);
        end

%--------------------------------------------------------------------------
%       Clean up data and fix problems
%--------------------------------------------------------------------------

        function prepareFixes(obj)
            % prepareFixes(obj) sets up structures to find and fix problems
            arguments
                obj (1,:)
            end
            obj.checkForError('dataNotLoaded',{obj})
            for o=1:numel(obj)
                obj(o).proc.errorhandling = {};
                fc = obj(o).get('frameCount');
                if strcmp(obj(o).ctrl.tracking,'hat')
                    for i=1:obj(o).getLength
                        obj(o).data(i).bad1 = false(fc(i),1);
                        obj(o).data(i).bad2 = false(fc(i),2);
                    end
                elseif strcmp(obj(o).ctrl.tracking,'cog')
                    for i=1:obj(o).getLength
                        obj(o).data(i).bad1 = false(fc(i),1);
                    end
                end
                switch cell2mat(obj(o).getInfo('type'))
                    case 'SingleCaro'
                        obj(o).proc.protectedRois = {'disk'};
                    case 'DoubleCaro'
                        obj(o).proc.protectedRois = {'disk1','disk2'};
                end
            end
        end

        function fixOutOfArena(obj,params)
            % fixOutOfArena(obj,params) marks fly positions (head and tail)
            % when they are outside of the arena
            arguments
                obj (1,:)
                params.aggressiveness = 0
            end
            obj.checkForError('dataNotLoaded',{obj})
            m.name = 'fixOutOfArena';
            m.params = params;
            for o=1:numel(obj)
                center = obj(o).ctrl.rois.arena;
                if strcmp(obj(o).ctrl.tracking,'hat')
                    cvals = hypot(obj(o).get('xvals')-center.x,obj(o).get('yvals')-center.y);
                    idx = cvals > center.rmax-params.aggressiveness;
                    obj(o).assignDataInPortions('bad2',idx,obj(o).get('frameCount'));
                elseif strcmp(obj(o).ctrl.tracking,'cog')
                    cvals = obj(o).distanceToTargetCenter('arena');
                    idx = cvals > center.rmax-params.aggressiveness;
                    obj(o).assignDataInPortions('bad1',idx,obj(o).get('frameCount'));
                end
                obj(o).proc.errorhandling{end+1} = m;
            end
        end

        function findJumps(obj,targets)
            % findJumps(obj,targets) finds frames with movement distance
            % above a given threshold. 'targets' is a list of protected
            % ROIs
            arguments
                obj (1,:)
                targets {CaroInterface.mustBeRoiList(targets,obj)} = obj(1).proc.protectedRois
            end
            obj.checkForError('dataNotLoaded',{obj})
            m.name = 'findJumps';
            m.params = targets;
            for o=1:numel(obj)
                l = [0; obj(o).distanceMoved>obj(o).conf.maxMoveDist];
                b = any([obj(o).get('bad1') l],2);
                oat = false(obj(o).getDataLength,numel(targets));
                for t=1:numel(targets)
                    ot = obj(o).onTarget(targets{t});
                    oat(:,t) = and(ot,circshift(ot,1));
                end
                bc = and(b,~any(oat,2));
                obj(o).assignDataInPortions('bad1',bc,obj(o).fc);
                m.threshold = obj(o).conf.maxMoveDist;
                obj(o).proc.errorhandling{end+1} = m;
            end
        end

        function keepJumpsOnly(obj)
            % keepJumpsOnly(obj) keeps only points around jumping events
            % (for diagnostic reasons)
            arguments
                obj (1,1)
            end
            obj.checkForError('dataNotLoaded',{obj})
            d = obj.distanceMoved;
            l = d>obj.conf.maxMoveDist;
            idx = find(l);
            idx2 = idx-1;
            idx_all = unique([idx idx2]);
            x = obj.get('x');
            y = obj.get('y');
            xn = nan(numel(x),1);
            yn = nan(numel(x),1);
            xn(idx_all) = x(idx_all);
            yn(idx_all) = y(idx_all);
            fc = obj.get('frameCount');
            obj.assignDataInPortions('x',xn,fc);
            obj.assignDataInPortions('y',yn,fc);
        end

        function calculateXY(obj)
            % calculateXY(obj) calculates x and y values
            arguments
                obj (1,:)
            end
            obj.checkForError('dataNotLoaded',{obj})
            m.name = 'calculateXY';
            for o=1:numel(obj)
                if strcmp(obj(o).ctrl.tracking,'cog')
                    fprintf('skipping %s for cog dataset\n',m.name)
                    continue
                end
                fc = obj(o).get('frameCount');
                idx = obj(o).get('bad2');
                xv = obj(o).get('xvals');
                yv = obj(o).get('yvals');
                xv(idx) = NaN;
                yv(idx) = NaN;
                obj(o).assignDataInPortions('x',mean(xv,2,'omitnan'),fc);
                obj(o).assignDataInPortions('y',mean(yv,2,'omitnan'),fc);
                obj(o).proc.errorhandling{end+1} = m;
            end
        end

        function applyBad(obj)
            % applyBad(obj) sets problematic points to NaN, to be fixed
            % later
            arguments
                obj (1,:)
            end
            obj.checkForError('dataNotLoaded',{obj})
            m.name = 'applyBad';
            for o=1:numel(obj)
                x = obj(o).get('x');
                y = obj(o).get('y');
                idx = obj(o).get('bad1');
                x(idx) = NaN;
                y(idx) = NaN;
                obj(o).assignDataInPortions('x',x,obj(o).fc);
                obj(o).assignDataInPortions('y',y,obj(o).fc);
                obj(o).proc.errorhandling{end+1} = m;
            end
        end

        function finalizeFixes(obj)
            % finalizeFixes(obj) removes data not needed after fix
            arguments
                obj (1,:)
            end
            obj.checkForError('dataNotLoaded',{obj})
            m.name = 'finalizeFixes';
            for o=1:numel(obj)
                obj(o).removeField('bad1');
                if strcmp(obj(o).ctrl.tracking,'hat')
                    obj(o).removeField('bad2');
                    obj(o).removeField('xvals');
                    obj(o).removeField('yvals');
                    obj(o).removeField('qvals');
                end
                obj(o).proc.errorhandling{end+1} = m;
            end
        end

        function runFixes(obj,params)
            % runFixes(obj,params) runs the standard set of repair options
            arguments
                obj (1,:)
                params.aggressiveness = 0;
            end
            obj.checkForError('dataNotLoaded',{obj})
            obj.prepareFixes;
            obj.fixOutOfArena('aggressiveness',params.aggressiveness);
            obj.calculateXY;
            obj.findJumps;
            obj.applyBad;
            obj.finalizeFixes;
        end

    end

    methods %override

        function v = get(obj,varargin)
            if isempty(obj.data)
                disp('This list does not contain data!')
                v = [];
            else
                if nargin > 1
                    v = get@BasicDataHandleList(obj,varargin{:});
                else
                    get@BasicDataHandleList(obj);
                end
            end
        end

    end

    methods (Access = protected)

        function checkForError(obj,errorname,params)
            msgID = '';
            switch errorname
                case 'mixedExperimentTypes'
                    expTypes = unique(obj.getInfo('type'));
                    if numel(expTypes) ~= params{1}
                        msgID = 'CaroDataLoader:mixedExperimentTypes';
                        msgTXT = 'You can not process different experiment types at the same time.';
                    end
                case 'dataNotLoaded'
                    if ~all(arrayfun(@(o) o.proc.dataLoaded,params{1}))
                        msgID = 'CaroDataLoader:dataNotLoaded';
                        msgTXT = 'Data has not been loaded yet. Please run loadData first';
                    end
                otherwise
                    checkForError@BasicDataHandleList(obj,errorname,params)
            end
            if ~isempty(msgID)
                throwAsCaller(MException(msgID,[msgID '. ' msgTXT]))
            end
        end

    end

    methods (Static)

        function [cdl,errormsg] = LoadMetaDataFromFolder(datapath)
            % [cdl,errormsg] = LoadMetaDataFromFolder(datapath) loads only
            % the metadata from a folder
            arguments
                datapath char {mustBeFolder}
            end
            cdl = CaroDataLoader([]);
            fparts = strsplit(datapath,filesep);
            fullname = fparts{end};
            cdl.info.fullName = fullname;
            cdl.info.dataPath = datapath;
            cdl.info.metaDataFile = '';
            errormsg = [];
            try
                cdl.info.mp4Count = numel(dir([datapath '\*.mp4']));
                cdl.info.h264Count = numel(dir([datapath '\*.h264']));
                cdl.info.h5Count = numel(dir([datapath '\*.h5']));
                cdl.info.csvCount = numel(dir([datapath '\*.csv']));
                csvname = sprintf('%s\\%s.csv',datapath,fullname);
                jsonname = sprintf('%s\\experiment_settings.json',datapath);
                if exist(csvname,'file')
                    cdl.info.metaDataFile = csvname;
                    cdl.loadMetaData(csvname);
                    cdl.info.csvCount = cdl.info.csvCount-1;
                elseif exist(jsonname,'file')
                    cdl.info.metaDataFile = jsonname;
                    cdl.loadMetaData(jsonname);
                else
                    errormsg = ['No metadata found in ',fullname];
                    disp(errormsg)
                end
            catch ME
                errormsg = [ME.identifier ' in ',fullname];
            end
            cdl.proc.loading{end+1} = 'LoadFromFolder';
        end

        function cdl = LoadExperimentFromFolder(datapath)
            % cdl = LoadExperimentFromFolder(datapath) loads an experiment
            % from a folder
            arguments
                datapath char {mustBeFolder}
            end
            cdl = CaroDataLoader.LoadMetaDataFromFolder(datapath);
            if isempty(cdl.info.metaDataFile)
                cdl = [];
                return
            end
            cdl.loadData
        end

        function [cdl,errorlog] = ExtractMetaDataFromBaseFolder(datapath)
            % [cdl,errorlog] = ExtractMetaDataFromBaseFolder(datapath)
            % extracts metadata
            arguments
                datapath char {mustBeFolder}
            end
            disp('Loading metadata, please wait...')
            folders = dir(datapath);
            folders = folders([folders.isdir]);
            folders = folders(3:end);
            cdl = CaroDataLoader.empty(numel(folders),0);
            errorlog = {};
            for i=1:numel(folders)
                [cdl(i),errormsg] = CaroDataLoader.LoadMetaDataFromFolder(fullfile(datapath,folders(i).name));
                cdl(i).proc.loading{end+1} = 'ExtractMetaDataFromBaseFolder';
                if ~isempty(errormsg)
                    errorlog(end+1) = {errormsg}; %#ok<AGROW>
                end
            end
            %remove empty/faulty experiments
            if numel(cdl) > 0
                cdl(cellfun(@isempty,cdl.getInfo('metaDataFile'))) = [];
            else
                warning('No experiments found!')
            end
        end

        function cdhl = LoadExperimentQuick(datapath)
            % cdhl = LoadExperimentQuick(datapath) loads a single
            % experiment from 'datapath' and returns a 
            % (Double)CaroDataHandleList after applying the standard repair
            % options
            arguments
                datapath char {mustBeFolder}
            end
            cdl = CaroDataLoader.LoadExperimentFromFolder(datapath);
            if isempty(cdl)
                cdhl = [];
                return
            end
            cdl.runFixes;
            cdhl = cdl.toCaroDataHandleList;
        end

    end

end