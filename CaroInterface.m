classdef CaroInterface < handle
%CAROINTERFACE Helper class for CaroDataHandleList and CaroDataLoader
%   Defines basic concepts like ROIs
%
%   Author: Tilman Triphan, tilman.triphan@uni-leipzig.de
%   License: GNU General Public License v3.0
%   Copyright (C) 2017-2023  Tilman Triphan
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
        ctrl %physical parameters (frame rate, image resolution, ROI positions)
        conf %configuration, settings or preferences
    end

    methods

        function obj = CaroInterface()
            obj.setDefaultValues()
        end

        function n = getDataLength(obj,indices)
            % n = getDataLength(obj,indices) returns the total number of
            % frames or for a subset if parameter indices is set
            arguments
                obj (1,1)
                indices = ':'
            end
            n = numel(obj.get('x',indices));
        end

        function d = fc(obj,indices)
            % d = fc(obj,indices) returns the number of frames for each
            % video or a subset if parameter indices is set
            arguments
                obj (1,1)
                indices {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
            end
            d = obj.get('frameCount',indices);
        end

        function v = getRoiList(obj)
            % v = getRoiList(obj) returns a list of all available ROIs
            v = obj.conf.roilist;
        end

        function v = getRoiCode(obj,target)
            % v = getRoiCode(obj,target) returns a list and order of all active ROIs
            if nargin < 2
                v = obj.conf.roicode;
            else
                v = find(ismember(obj.conf.roilist,target));
            end
        end

        function v = getActiveRois(obj)
            % v = getActiveRois(obj) returns a list of all active ROIs
            arguments
                obj (1,1)
            end
            v = obj.conf.roilist(obj.conf.roicode);
        end

        function addCustomInfoGroup(obj,fieldname,sourcename,groups)
            % addCustomInfoGroup(obj,fieldname,sourcename,groups) adds
            % fields to the info struct to put experiments into groups
            arguments
                obj (1,:)
                fieldname (1,:) char
                sourcename (1,:) char
                groups (1,:)
            end
            if isfield(obj.getInfo,sourcename)
                v = obj.getInfo(sourcename);
                v = [v{:}];
            else
                v = obj.(sourcename);
            end
            d = nan(numel(obj),1);
            for i=1:numel(groups)
                idx = v>=groups{i}(1) & v<groups{i}(end);
                d(idx) = i;
            end
            obj.addInfo(fieldname,num2cell(d));
        end

        function l = isNaN(obj,indices)
            % l = isNaN(obj,indices) returns a logical array with true
            % values indicating NaNs for x (i.e. failures in fly detection)
            arguments
                obj (1,1)
                indices {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
            end
            l = isnan(obj.get('x',indices));
        end

        function [theta, rho] = toPolar(obj,target,indices)
            % [theta, rho] = toPolar(obj,target,indices) returns polar
            % coordinates in respect to the given target ROI
            arguments
                obj (1,1)
                target {CaroInterface.mustBeRoiName(target,obj)} = 'arena'
                indices {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
            end
            xloc = obj.get('x',indices) - obj.ctrl.rois.(target).x;
            yloc = obj.get('y',indices) - obj.ctrl.rois.(target).y;
            [theta, rho] = cart2pol(xloc,yloc);
        end

        function d = distanceMoved(obj,indices)
            % d = distanceMoved(obj,indices) returns the distance moved
            % (active or passive) between two frames
            arguments
                obj (1,1)
                indices {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
            end
            d = [hypot(diff(obj.get('x',indices)),diff(obj.get('y',indices))) NaN];
        end

        function d = distanceWalked(obj,indices)
            % d = distanceWalked(obj,indices) returns the distance walked
            % actively by the fly between two frames. Movement on the
            % carousel is masked by NaNs
            arguments
                obj (1,1)
                indices {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
            end
            d = obj.distanceMoved(indices);
            d(obj.onTarget('disk',indices)) = NaN;
        end

        function d = sumDistanceMoved(obj,indices)
            % d = sumDistanceMoved(obj,indices) returns the sum of the
            % distance moved (active or passive) for each experiment in the
            % list
            arguments
                obj (1,:)
                indices {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
            end
            d = arrayfun(@(o) sum(o.distanceMoved(indices),'omitnan'),obj)';
        end

        function d = sumDistanceWalked(obj,indices)
            % d = sumDistanceWalked(obj,indices) returns the sum of the
            % distance walked actively be the fly for each experiment in
            % the list
            arguments
                obj (1,:)
                indices {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
            end
            d = arrayfun(@(o) sum(o.distanceWalked(indices),'omitnan'),obj)';
        end

        function d = meanDistanceMoved(obj,indices)
            % d = meanDistanceMoved(obj,indices) returns the mean of
            % distanceMoved for each experiment in the list
            arguments
                obj (1,:)
                indices {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
            end
            d = arrayfun(@(o) mean(o.distanceWalked(indices),'omitnan'),obj)';
        end

        function d = meanDistanceWalked(obj,indices)
            % d = meanDistanceWalked(obj,indices) returns the mean of
            % distanceWalked for each experiment in the list
            arguments
                obj (1,:)
                indices {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
            end
            d = arrayfun(@(o) mean(o.distanceWalked(indices),'omitnan'),obj)';
        end

        function d = distanceToTargetCenter(obj,target,indices)
            % d = distanceToTargetCenter(obj,target,indices) returns the
            % distance to the center of the given target ROI
            arguments
                obj (1,1)
                target {CaroInterface.mustBeRoiName(target,obj)} = 'arena'
                indices {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
            end
            [~,d] = obj.toPolar(target,indices);
        end

        function l = isActive(obj,indices,params)
            % l = isActive(obj,indices,params) returns a logical array with
            % 'true' indicating all frames where the fly moved a larger
            % distance than minActivity
            arguments
                obj (1,1)
                indices {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
                params.minActivity (1,1) double = obj.conf.minActivity;
            end
            d = obj.distanceWalked(indices);
            l = d > params.minActivity;
        end

        function l = onTarget(obj,target,indices)
            % l = onTarget(obj,target,indices) returns a logical array with
            % 'true' for all frames where the fly is in the 'target' ROI
            arguments
                obj (1,1)
                target {CaroInterface.mustBeRoiName(target,obj)}
                indices {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
            end
            t = obj.ctrl.rois.(target);
            switch (t.type)
                case 'circ'
                    d = obj.distanceToTargetCenter(target,indices);
                    l = d >= t.rmin & d <= t.rmax;
                case 'func'
                    l = obj.(t.func)(indices);
                case 'comb'
                    inc = t.inc;
                    l0 = false(obj.getDataLength(indices),numel(inc));
                    for i=1:numel(inc)
                        l0(:,i) = obj.onTarget(inc{i},indices);
                    end
                    switch (t.rule)
                        case {'and','all'}
                            l = all(l0,2);
                        case {'or','any'}
                            l = any(l0,2);
                        case {'not'}
                            l0(:,2:end) = ~l0(:,2:end);
                            l = all(l0,2);
                    end
            end
        end

        function l = isPointOnTarget(obj,x,y,target)
            % l = isPointOnTarget(obj,x,y,target) returns 'true' if a point
            % with the coordinates 'x' and 'y' is in the 'target' ROI
            arguments
                obj (1,1)
                x (1,1) single {CaroInterface.mustBeValidXPosition(x,obj)}
                y (1,1) single {CaroInterface.mustBeValidYPosition(y,obj)}
                target (1,:) char {CaroInterface.mustBeRoiName(target,obj)}
            end
            roi = obj.ctrl.rois.(target);
            if ~strcmp(roi.type,'circ')
                error('CaroInterface:incorrectRoiType','isPointOnTarget only handles rois of type circ')
            end
            xloc = x - roi.x;
            yloc = y - roi.y;
            [~, rho] = cart2pol(xloc,yloc);
            l = (rho >= roi.rmin) && (rho <= roi.rmax);
        end

        function r = getRoiForPoint(obj,x,y)
            % r = getRoiForPoint(obj,x,y) returns the respective ROI (or
            % ROIs) for a set of x and y coordinates
            arguments
                obj (1,1)
                x (1,1) single {CaroInterface.mustBeValidXPosition(x,obj)}
                y (1,1) single {CaroInterface.mustBeValidYPosition(y,obj)}
            end
            rois = obj.ActiveRois;
            r = {};
            for f = 1:numel(rois)
                if obj.isPointOnTarget(x,y,rois{f})
                    r{end+1} = rois{f}; %#ok<AGROW>
                end
            end
        end

        function f = fractionOnTarget(obj,target,indices)
            % f = fractionOnTarget(obj,target,indices) returns the fraction
            % of frames where the fly was on the 'target' ROI
            arguments
                obj (1,:)
                target {CaroInterface.mustBeRoiName(target,obj)}
                indices {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
            end
            f = arrayfun(@(o) mean(o.onTarget(target,indices),'omitnan'),obj);
        end

        function f = fractionNaN(obj,indices)
            % f = fractionNaN(obj,indices) returns the fraction of frames
            % where the fly was not detected
            arguments
                obj (1,:)
                indices {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
            end
            f = arrayfun(@(o) mean(o.isNaN(indices),'omitnan'),obj);
        end

        function f = fractionActive(obj,indices)
            % f = fractionActive(obj,indices) returns the fraction of
            % frames where the fly is active, i.e. moving above a given
            % threshold
            arguments
                obj (1,:)
                indices {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
            end
            f = arrayfun(@(o) mean(o.isActive(indices),'omitnan'),obj);
        end

        function invalidateData(obj,frames)
            % invalidateData(obj,frames) sets the x values given in
            % 'frames' to NaN
            arguments
                obj (1,1)
                frames {CaroInterface.mustBeValidFrameIndices(frames,obj)}
            end
            x = obj.get('x');
            x(frames) = NaN;
            obj.assignDataInPortions('x',x,obj.fc);
        end

        function invalidateRegion(obj,x,y,r)
            % invalidateRegion(obj,x,y,r) invalidates (i.e. sets to NaN)
            % all points in a circle defined by the center coordinates 'x'
            % and 'y' and the radius 'r'
            arguments
                obj (1,1)
                x (1,1) single {CaroInterface.mustBeValidXPosition(x,obj)}
                y (1,1) single {CaroInterface.mustBeValidYPosition(y,obj)}
                r (1,1) double
            end
            obj.addRoi('eraser','x',x,'y',y,'rmax',r,'type','circ');
            obj.invalidateData(obj.onTarget('eraser'));
            obj.removeRoi('eraser');
        end


        function setXY(obj,x,y,frames)
            % setXY(obj,x,y,frames) sets the x and y values for the given
            % frame numbers
            arguments
                obj (1,1)
                x (1,1) single {CaroInterface.mustBeValidXPosition(x,obj)}
                y (1,1) single {CaroInterface.mustBeValidYPosition(y,obj)}
                frames {CaroInterface.mustBeValidFrameIndices(frames,obj)}
            end
            xv = obj.get('x');
            yv = obj.get('y');
            xv(frames) = x;
            yv(frames) = y;
            obj.assignDataInPortions('x',xv,obj.fc);
            obj.assignDataInPortions('y',yv,obj.fc);
        end

        function setDefaultValues(obj,filename)
            % setDefaultValues(obj,filename) reads default parameters from
            % a json file and sets them
            arguments
                obj (1,1)
                filename {mustBeFile} = fullfile(fileparts(which('/CaroInterface.m')),'defaultSettings.json')
            end
            obj.conf = jsondecode(fileread(filename));
        end

        function updateRoiConf(obj)
            % updateRoiConf(obj) updates ROI handling for plotting
            arguments
                obj (1,1)
            end
            an = fieldnames(obj.ctrl.rois);
            ac = [];
            as = [];
            ap = [];
            obj.conf.roilist = an;
            obj.conf.cmapdefault = nan(numel(an),3);
            for i=1:numel(an)
                obj.conf.cmapdefault(i,:) = obj.ctrl.rois.(an{i}).color;
                if obj.ctrl.rois.(an{i}).use
                    ac = [ac i]; %#ok<AGROW>
                    as = [as obj.ctrl.rois.(an{i}).sort]; %#ok<AGROW>
                    ap = [ap obj.ctrl.rois.(an{i}).prio]; %#ok<AGROW>
                end
            end
            [~,idx] = sort(as);
            obj.conf.roicode = ac(idx);
            [~,idx] = sort(ap);
            obj.conf.areaprio = ac(idx);
        end

        function r = addRoi(obj,name,params)
            % r = addRoi(obj,name,params) adds a new ROI or updates a
            % existing ROI
            arguments
                obj (1,1)
                name (1,:) char
                params.x (1,1) double = -1
                params.y (1,1) double = -1
                params.rmin (1,1) double = 0
                params.rmax (1,1) double = -1
                params.color (1,3) double = [0 0 0]
                params.comment (1,:) char = ''
                params.func (1,:) char = '' %ismethod
                params.prio (1,1) double = 0
                params.show (1,1) logical = false
                params.sort (1,1) double = 0
                params.type (1,:) char {mustBeMember(params.type,{'circ','func','comb','spec'})}
                params.use (1,1) logical = false
            end

            r.type = params.type;
            r.color = params.color;
            r.prio = params.prio;
            r.sort = params.sort;
            switch params.type
                case 'circ'
                    r.x = params.x;
                    r.y = params.y;
                    r.rmin = params.rmin;
                    r.rmax = params.rmax;
                case 'func'
                    if ~ismethod(obj,params.func)
                        disp('Not a valid method')
                    else
                        r.func = params.func;
                    end
            end
            r.comment = params.comment;
            r.use = params.use;
            r = orderfields(r);
            obj.ctrl.rois.(name) = r;
            obj.updateRoiConf;
        end

        function removeRoi(obj,target)
            % removeRoi(obj,target) removes a ROI
            arguments
                obj (1,1)
                target (1,:) char {CaroInterface.mustBeRoiName(target,obj)}
            end
            obj.ctrl.rois = rmfield(obj.ctrl.rois,target);
            obj.updateRoiConf;
        end

        function assignDataInPortions(obj,fieldname,data,portions,indices)
            % assignDataInPortions assigns data to a field into diffently
            % sized portions (e.g. for updating videos with variable frame
            % numbers). The number of portions must be equal to the number
            % of assignees and the total size must be valid as well.
            arguments
                obj (1,1)
                fieldname {CaroInterface.mustBeValidDataField(fieldname,obj)}
                data
                portions
                indices {CaroInterface.mustBeValidIndex(indices,obj)} = 1:obj.getLength
            end
            idx = find(indices);
            if size(data,1) ~= sum(portions)
                error('The total sum of all portions does not match the data size')
            end
            if numel(idx) ~= numel(portions)
                error('The number of portions does not match the number of assignees')
            end
            offset = 0;
            for i=1:numel(idx)
                d = obj.getData(idx(i));
                d.(fieldname) = data(offset+1:offset+portions(i),:);
                obj.setData(d,i);
                offset = offset + portions(i);
            end
        end


    end

    methods (Static, Hidden = true)

        function mustBeRoiName(target,obj)
            % Validate that 'target' is a valid ROI name.
            if ~all(arrayfun(@(o) ismember(target,o.conf.roilist),obj))
                eid = "mustBeRoiName:invalidRoiName";
                msg = "'target' must be a valid ROI name for all items in list.";
                throwAsCaller(MException(eid,msg))
            end
        end

        function mustBeRoiList(targets,obj)
            % Validate that 'targets' only contains valid ROI names.
            if ~all(cellfun(@(z) all(arrayfun(@(o) ismember(targets,o.conf.roilist),z)),obj))
                eid = "mustBeRoiList:invalidRoiName";
                msg = "'targets' must only contain valid ROI names.";
                throwAsCaller(MException(eid,msg))
            end
        end

        function mustBeListsInCell(groups)
            % Validate that 'groups' is a cell array containing only
            % CaroDataHandleLists
            if ~all(cellfun(@(x)isa(x,'CaroDataHandleList'),groups))
                eid = "mustBeListsInCell:invalidCellContents";
                msg = "'groups' must only contain CaroDataHandleLists.";
                throwAsCaller(MException(eid,msg))
            end
        end

        function mustBeValidIndex(indices,obj)
            % Validate that 'indices' are valid indices for the data.
            % 'indices' can be either logical or numerical. The valid range
            % depends on the .conf.indexMode (either 'frame' or 'video').
            if strcmp(indices,':')
                return
            elseif islogical(indices)
                if strcmp(obj(1).conf.indexMode,'frame') && (numel(indices) == obj(1).getDataLength)
                    return
                elseif strcmp(obj(1).conf.indexMode,'video') && (numel(indices) == obj(1).getLength)
                    return
                else
                    eid = "mustBeValidIndex:invalidLogicalIndex";
                    msg = "'indices' must be valid indices.";
                    throwAsCaller(MException(eid,msg))
                end
            elseif all(indices > 0) && all(indices == floor(indices))
                if strcmp(obj(1).conf.indexMode,'frame') && all(indices <= obj(1).getDataLength)
                    return
                elseif strcmp(obj(1).conf.indexMode,'video') && all(indices <= obj(1).getLength)
                    return
                else
                    eid = "mustBeValidIndex:invalidNumericalIndex";
                    msg = "'indices' must be valid indices.";
                    throwAsCaller(MException(eid,msg))
                end
            else
                eid = "mustBeValidIndex:invalidIndex";
                msg = "'indices' must be valid indices.";
                throwAsCaller(MException(eid,msg))
            end
        end

        function mustBeValidFrameIndices(frames,obj)
            % Validate that 'frames' are valid frame indices for the data.
            % 'frames' can be either logical or numerical.
            if islogical(frames) && (numel(frames) == obj.getDataLength)
                return
            elseif all(frames > 0) && all(frames <= obj.getDataLength)
                return
            else
                eid = "mustBeValidFrames:invalidFrameIndex";
                msg = "'frames' must be valid frame indices.";
                throwAsCaller(MException(eid,msg))
            end
        end

        function mustBeValidDataField(fieldName,obj)
            % Validate that 'fieldName' is a valid fieldname of the data
            % structure
            if ~isfield(obj.getData,fieldName)
                eid = "mustBeValidDataField:invalidFieldName";
                msg = "'fieldName' must be a data fieldname.";
                throwAsCaller(MException(eid,msg))
            end
        end

        function mustBeValidInfoField(fieldName,obj)
            % Validate that 'fieldName' is a valid fieldname of the info
            % structure
            if ~isfield(obj.getInfo,fieldName)
                eid = "mustBeValidInfoField:invalidFieldName";
                msg = "'fieldName' must be a info fieldname.";
                throwAsCaller(MException(eid,msg))
            end
        end

        function mustBeValidInfoFieldOrMethod(fieldName,obj)
            % Validate that 'fieldName' is a valid fieldname of the info
            % structure
            if ~isfield(obj.getInfo,fieldName) && ~ismethod(obj,fieldName)
                eid = "mustBeValidInfoFieldOrMethod:invalidFieldName";
                msg = "'fieldName' must be a info fieldname or method.";
                throwAsCaller(MException(eid,msg))
            end
        end

        function mustBeValidDataFieldOrMethod(fieldName,obj)
            % Validate that 'fieldName' is a valid fieldname of the data
            % structure or a valid method
            if ~any([isfield(obj.getData,fieldName),ismethod(obj,fieldName)])
                eid = "mustBeValidDataFieldOrMethod:invalidFieldOrMethodName";
                msg = "'fieldName' must be a data fieldname or method.";
                throwAsCaller(MException(eid,msg))
            end
        end

        function mustBeIdenticalInfoValue(fieldName,obj)
            % Validate that the value in 'fieldName' is identical for all
            % elements in the list
            if ~numel(unique(obj.getInfo(fieldName))) == 1
                eid = "mustBeIdenticalInfoValue:inconsistentValue";
                msg = ['the value in ',fieldName,' must be identical for all elements in the list.'];
                throwAsCaller(MException(eid,msg))
            end
        end

        function mustBeValidXPosition(x,obj)
            % Validate that the value in 'x' is a valid x position between
            % 1 and imageResX
            if x >= 1 && x <= obj.ctrl.imageResX
                return
            else
                eid = "mustBeValidXPosition:invalidXPosition";
                msg = "'x' must be a valid x position.";
                throwAsCaller(MException(eid,msg))
            end
        end

        function mustBeValidYPosition(y,obj)
            % Validate that the value in 'y' is a valid y position between
            % 1 and imageResY
            if y >= 1 && y <= obj.ctrl.imageResY
                return
            else
                eid = "mustBeValidYPosition:invalidYPosition";
                msg = "'y' must be a valid y position.";
                throwAsCaller(MException(eid,msg))
            end
        end

    end

end