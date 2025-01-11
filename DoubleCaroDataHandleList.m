classdef DoubleCaroDataHandleList < CaroDataHandleList
%DOUBLECARODATAHANDLELIST  Class to handle DrosoCarousel experiments with
%   two carousels
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

    methods

        function obj = DoubleCaroDataHandleList(data)
            if nargin > 0
                obj.data = data;
            end
        end

        function l = isDisk1Moving(obj,indices)
            % l = isDisk1Moving(obj,indices) returns true if the first
            % carousel is turning
            arguments
                obj (1,1)
                indices {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
            end
            if any(contains(obj.proc.loading,'loadDisk'))
                l = obj.get('rot1',indices);
            else
                disp('run loadDisk first!')
            end
        end

        function l = isDisk2Moving(obj,indices)
            % l = isDisk2Moving(obj,indices) returns true if the second
            % carousel is turning
            arguments
                obj (1,1)
                indices {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
            end
            if any(contains(obj.proc.loading,'loadDisk'))
                l = obj.get('rot2',indices);
            else
                disp('run loadDisk first!')
            end
        end

        function setDiskRules(obj,rules)
            % setDiskRules(obj,rules) allows switching the analysis rules. 
            % 'diskOff' uses the ROIs 'disk1' and 'disk2' and only looks if
            % the fly is sitting on one or the other carousel but does not
            % care if the carousel is turning.
            % 'diskOn' uses the ROIs 'diskturn1' and 'diskturn2' as well as
            % 'diskstop1' and 'diskstop2'.
            arguments
                obj (1,1)
                rules (1,:) char {mustBeMember(rules,{'diskOn','diskOff'})} = 'diskOff'
            end
            switch rules
                case 'diskOn'
                    obj.ctrl.rois.turn1.use = false;
                    obj.ctrl.rois.turn2.use = false;
                    obj.ctrl.rois.diskturn1.use = true;
                    obj.ctrl.rois.diskturn2.use = true;
                    obj.ctrl.rois.diskstop1.use = true;
                    obj.ctrl.rois.diskstop2.use = true;
                    obj.ctrl.rois.diskstop1.prio = -50;
                    obj.ctrl.rois.diskstop2.prio = -50;
                    obj.ctrl.rois.disk1.use = false;
                    obj.ctrl.rois.disk2.use = false;
                case 'diskOff'
                    obj.ctrl.rois.turn1.use = false;
                    obj.ctrl.rois.turn2.use = false;
                    obj.ctrl.rois.diskturn1.use = false;
                    obj.ctrl.rois.diskturn2.use = false;
                    obj.ctrl.rois.diskstop1.use = false;
                    obj.ctrl.rois.diskstop2.use = false;
                    obj.ctrl.rois.disk1.use = true;
                    obj.ctrl.rois.disk2.use = true;
            end
            obj.proc.rules = rules;
            obj.updateRoiConf();
        end

        %TODO: remove
        function h = plotDiskTransitionCharts(obj,indices,params)
            % h = plotDiskTransitionCharts(obj,indices,params) plots pie
            % charts for the transitions between left and right carousel
            % starting and stopping
            arguments
                obj (1,:)
                indices = ':'
                params.normalize (1,1) logical = true
            end
            if ~all(arrayfun(@(o) ismember('loadDiskData',o.proc.loading),obj))
                disp('load disk turn info first')
                h = [];
                return
            end
            h = figure;
            al = obj(1).conf.roilist;
            rl = {'diskturn1','diskstop1','diskturn2','diskstop2'};
            colormap(obj(1).conf.cmapdefault)
            hc = nan(numel(al),numel(obj),4);

            for o=1:numel(obj)
                e = obj(o).createEthogram(indices);
                e(e==0) = find(ismember(al,'arena'));
                for i=1:numel(rl)
                    v = e((obj(o).findTransitions(rl{i},indices)));
                    if strcmp(al(e(1)),rl{i})
                        v = v(1:2:end);
                    else
                        v = v(2:2:end);
                    end
                    n = nan(numel(al),1);
                    for j=1:numel(al)
                        n(j) = sum(v==j);
                    end
                    if params.normalize
                        hc(:,o,i) = n/sum(n);
                    else
                        hc(:,o,i) = n;
                    end
                end
            end
            for i=1:numel(rl)
                labels = al;
                labels([3 5 8 10:14 19 20]) = {''};
                subplot(2,2,i)
                pie(mean(hc(:,:,i),2,'omitnan'),labels)
                title(rl{i})
            end
        end

        function [h,t] = plotDiskTransitionDifference(obj,indices,params)
            % [h,t] = plotDiskTransitionDifference(obj,indices,params)
            arguments
                obj (1,:)
                indices = ':'
                params.excludeNan (1,1) logical = true
            end
            if ~all(arrayfun(@(o) ismember('loadDiskData',o.proc.loading),obj))
                disp('load disk turn info first')
                h = [];
                return
            end
            h = figure;
            dt1v = nan(numel(obj),1);
            ds1v = nan(numel(obj),1);
            dt2v = nan(numel(obj),1);
            ds2v = nan(numel(obj),1);
            for o=1:numel(obj)
                al = obj(o).conf.roilist;
                e = obj(o).createEthogram(indices);
                e(e==0) = find(ismember(al,'arena'));

                dt1 = e((obj(o).findTransitions('diskturn1',indices)));
                if strcmp(al(e(1)),'diskturn1')
                    dt1 = dt1(1:2:end);
                else
                    dt1 = dt1(2:2:end);
                end
                if params.excludeNan
                    dt1v(o) = sum(dt1==find(ismember(al,'diskstop1')))/(numel(dt1)-sum(dt1==find(ismember(al,'nan'))));
                else
                    dt1v(o) = sum(dt1==find(ismember(al,'diskstop1')))/numel(dt1);
                end

                ds1 = e((obj(o).findTransitions('diskstop1',indices)));
                if strcmp(al(e(1)),'diskstop1')
                    ds1 = ds1(1:2:end);
                else
                    ds1 = ds1(2:2:end);
                end
                if params.excludeNan
                    ds1v(o) = sum(ds1==find(ismember(al,'diskturn1')))/(numel(ds1)-sum(ds1==find(ismember(al,'nan'))));
                else
                    ds1v(o) = sum(ds1==find(ismember(al,'diskturn1')))/numel(ds1);
                end

                dt2 = e((obj(o).findTransitions('diskturn2',indices)));
                if strcmp(al(e(1)),'diskturn2')
                    dt2 = dt2(1:2:end);
                else
                    dt2 = dt2(2:2:end);
                end
                if params.excludeNan
                    dt2v(o) = sum(dt2==find(ismember(al,'diskstop2')))/(numel(dt2)-sum(dt2==find(ismember(al,'nan'))));
                else
                    dt2v(o) = sum(dt2==find(ismember(al,'diskstop2')))/numel(dt2);
                end

                ds2 = e((obj(o).findTransitions('diskstop2',indices)));
                if strcmp(al(e(1)),'diskstop2')
                    ds2 = ds2(1:2:end);
                else
                    ds2 = ds2(2:2:end);
                end
                if params.excludeNan
                    ds2v(o) = sum(ds2==find(ismember(al,'diskturn2')))/(numel(ds2)-sum(ds2==find(ismember(al,'nan'))));
                else
                    ds2v(o) = sum(ds2==find(ismember(al,'diskturn2')))/numel(ds2);
                end
            end
            subplot(2,1,1)
            boxplot([dt1v ds1v])
            ylim([0 1])
            xticklabels({'diskstop1','diskturn1'})
            hold on
            for o=1:numel(obj)
                line([1 2],[dt1v(o) ds1v(o)],'Marker', '.', 'MarkerSize', 10)
            end
            title('disk1')
            subplot(2,1,2)
            boxplot([dt2v ds2v])
            ylim([0 1])
            xticklabels({'diskstop2','diskturn2'})
            hold on
            for o=1:numel(obj)
                line([1 2],[dt2v(o) ds2v(o)],'Marker', '.', 'MarkerSize', 10)
            end
            title('disk2')
            d = nan(numel(obj),4);
            d(1:numel(obj),1) = dt1v;
            d(1:numel(obj),2) = ds1v;
            d(1:numel(obj),3) = dt2v;
            d(1:numel(obj),4) = ds2v;
            t = array2table(d,'VariableNames',{'diskstop1','diskturn1','diskstop2','diskturn2'});
            try
                writetable(t,[pwd '\diskTransitionDifference.xlsx']);
            catch
                disp('failed to write table!')
            end
        end

        function [h,t] = plotDiskTransitionComparison(obj,cdhl,labels)
            % [h,t] = plotDiskTransitionComparison(obj,cdhl,labels)
            arguments
                obj (1,:)
                cdhl (1,:) DoubleCaroDataHandleList
                labels (1,2) cell = {'',''}
            end
            dt1_dt2_1 = obj.findSpecialTransitions('diskturn1','sortby','diskturn2','nofigure',true);
            dt2_dt1_1 = obj.findSpecialTransitions('diskturn2','sortby','diskturn1','nofigure',true);
            dt1_dt2_2 = cdhl.findSpecialTransitions('diskturn1','sortby','diskturn2','nofigure',true);
            dt2_dt1_2 = cdhl.findSpecialTransitions('diskturn2','sortby','diskturn1','nofigure',true);
            dt_1 = [dt1_dt2_1.count] + [dt2_dt1_1.count];
            dt_2 = [dt1_dt2_2.count] + [dt2_dt1_2.count];
            ds1_ds2_1 = obj.findSpecialTransitions('diskstop1','sortby','diskstop2','nofigure',true);
            ds2_ds1_1 = obj.findSpecialTransitions('diskstop2','sortby','diskstop1','nofigure',true);
            ds1_ds2_2 = cdhl.findSpecialTransitions('diskstop1','sortby','diskstop2','nofigure',true);
            ds2_ds1_2 = cdhl.findSpecialTransitions('diskstop2','sortby','diskstop1','nofigure',true);
            ds_1 = [ds1_ds2_1.count] + [ds2_ds1_1.count];
            ds_2 = [ds1_ds2_2.count] + [ds2_ds1_2.count];
            h = figure;
            vals = nan(max([numel(dt_1),numel(dt_2),numel(ds_1),numel(ds_2)]),4);
            vals(1:numel(dt_1),1) = dt_1';
            vals(1:numel(dt_2),2) = dt_2';
            vals(1:numel(ds_1),3) = ds_1';
            vals(1:numel(ds_2),4) = ds_2';
            boxplot(vals)
            t = array2table(vals);
            xticklabels({['on>on_',labels{1}],['on>on_',labels{2}],...
                ['off>off_',labels{1}],['off>off_',labels{2}]})
            t.Properties.VariableNames = xticklabels;
        end

        function findCarouselStays(obj,targets,params)
            % findCarouselStays(obj,targets,params) extracts events of the
            % fly moving between the turning carousels
            arguments
                obj (1,:)
                targets (1,:) cell {CaroInterface.mustBeRoiList(targets,obj)} = {'diskturn1','diskturn2'}
                params.threshold (1,1) double = 0.1
                params.maxgap (1,1) double = 2
                params.minstretch (1,1) double = 3
                params.postproc (1,1) logical = false
                params.allowonly (1,:) cell {CaroInterface.mustBeRoiList(params.allowonly,obj)} = {'ring1','ring2'}
                params.mintrans (1,1) double = 3
                params.plot (1,1) string {mustBeMember(params.plot,{'none','walkingTrace','positions'})} = 'walkingTrace'
                params.addframes (1,1) double = 100
            end
            for o=1:numel(obj)
                fc = obj(o).fc;
                lists = obj(o).split;
                tl = false(numel(lists),numel(targets));
                for t=1:numel(targets)
                    tl(:,t) = lists.fractionOnTarget(targets{t}) >= params.threshold;
                end
                im_backup = obj(o).conf.indexMode;
                obj(o).conf.indexMode = 'video';
                ot = all(tl,2);
                vids = find(ot);
                p = find(diff(vids)>params.maxgap);
                sp = [1; p+1];
                ep = [p; numel(vids)];
                idx = find(ep-sp+1 >= params.minstretch);
                for i=1:numel(idx)
                    vid_sp = vids(sp(idx(i)));
                    vid_ep = vids(ep(idx(i)));
                    fprintf('#:%d %d-%d %s\n',o,vid_sp,vid_ep,cell2mat(obj(o).getInfo('fullName')))
                    if params.postproc
                        % exclude transitions coming from the border roi
                        ap = obj(o).conf.roiprio;
                        obj(o).conf.roiprio = [17 18 4 6 7 9 2 15 16];
                        etho = obj(o).createEthogram(vid_sp:vid_ep+1);%get longer stretch
                        [~,~,ton1,~,~,~] = obj(o).getVisitMetrics('diskturn1',vid_sp:vid_ep+1);
                        [~,~,ton2,~,~,~] = obj(o).getVisitMetrics('diskturn2',vid_sp:vid_ep+1);
                        framecount = sum(obj(o).get('frameCount',(vid_sp:vid_ep)));
                        ton1(ton1>framecount) = [];
                        ton2(ton2>framecount) = [];
                        ton1(ton1==1) = [];
                        ton2(ton2==1) = [];
                        ac1 = etho(ton1-1);
                        ac2 = etho(ton2-1);
                        ac1(ac1==0) = 1;
                        ac2(ac2==0) = 1;
                        n1 = obj(o).conf.roilist(ac1)';
                        n2 = obj(o).conf.roilist(ac2)';
                        trevt = [];
                        trevt_nr = 0;
                        for ii=1:numel(ton1)
                            trevt_nr = trevt_nr +1;
                            trevt(trevt_nr).disk = 1;
                            trevt(trevt_nr).ton = ton1(ii);
                            trevt(trevt_nr).type = cell2mat(n1(ii));
                        end
                        for ii=1:numel(ton2)
                            trevt_nr = trevt_nr +1;
                            trevt(trevt_nr).disk = 2;
                            trevt(trevt_nr).ton = ton2(ii);
                            trevt(trevt_nr).type = cell2mat(n2(ii));
                        end
                        obj(o).conf.roiprio = ap;
                        [~,sort_idx] = sort([trevt.ton]);
                        trevt = trevt(sort_idx);
                        % find stretches of valid transitions
                        stretch_logical = ismember({trevt.type},params.allowonly);
                        disp(stretch_logical)
                        cc = bwconncomp(stretch_logical);
                        for c=1:numel(cc.PixelIdxList)
                            n = numel(cc.PixelIdxList{c});
                            if n >= params.mintrans
                                ids = cc.PixelIdxList{c};

                                trevt_run = trevt(ids);
                                global_frame_start = sum(fc(1:vid_sp-1));
                                for t=1:numel(trevt_run)
                                    plot_frame = trevt_run(t).ton-trevt_run(1).ton+params.addframes;
                                    global_frame = global_frame_start + plot_frame - params.addframes;
                                    fprintf('disk:%d frame:%d(l) %d(f) %d(g) from %s\n',trevt_run(t).disk,trevt_run(t).ton,plot_frame,global_frame,trevt_run(t).type)
                                end

                                if ~strcmp(params.plot,'none')
                                    cdhl = obj(o).createSubsetByIndex(vid_sp:vid_ep);
                                    cdhl.conf.indexMode = 'frame';
                                    indices = max(1,trevt(ids(1)).ton-params.addframes):min(trevt(ids(end)).ton+params.addframes,cdhl.getDataLength);
                                    switch params.plot
                                        case 'walkingTrace'
                                            cdhl.plotWalkingTrace(indices);
                                        case 'positions'
                                            cdhl.plotPositions('indices',indices);
                                    end
                                    view(0,30)
                                    h = gcf;
                                    gfs = global_frame_start + trevt_run(1).ton - params.addframes;
                                    gfe = global_frame_start + trevt_run(end).ton + params.addframes;
                                    h.Name = [h.Name sprintf(' vid:%d-%d frames:%d-%d',vid_sp,vid_ep,gfs,gfe)];
                                    fprintf('vid:%d:%d frames:%d:%d\n',vid_sp,vid_ep,gfs,gfe)
                                end
                            end
                        end
                    end
                    if ~strcmp(params.plot,'none') && ~ params.postproc
                        switch params.plot
                            case 'walkingTrace'
                                obj(o).plotWalkingTrace(vid_sp:vid_ep);
                            case 'positions'
                                obj(o).plotPositions('indices',vid_sp:vid_ep);
                        end
                        view(0,30)
                    end
                end
                obj(o).conf.indexMode = im_backup;
            end
        end

        function loadDiskData(obj,offs)
            % loadDiskData(obj,offs) loads information if the left or the
            % right carousel is turning
            arguments
                obj (1,1)
                offs (2,1) = [0.6,1.2]
            end
            try
                rot = cell(obj.getLength,2);
                for i=1:obj.getLength
                    fn = [cell2mat(obj.getInfo('dataPath')) '\' cell2mat(obj.get('fileID',i)) '.txt'];
                    t = readtable(fn);
                    d = t.Mean;
                    rot{i,1} = movmean(d(1:2:end) > offs(1),3) > 0.5;
                    rot{i,2} = movmean(d(2:2:end) > offs(2),3) > 0.5;
                end
                obj.addFieldAndValues('rot1',rot(:,1));
                obj.addFieldAndValues('rot2',rot(:,2));
                obj.proc.loading{end+1} = 'loadDiskData';
            catch
                disp('Failed to load disk data')
            end
        end
    end

    methods %override

        function d = distanceWalked(obj,indices)
            % d = distanceWalked(obj,indices) returns the distance walked
            % actively by the fly between two frames. Movement on the
            % active carousel is masked by NaNs
            arguments
                obj (1,1)
                indices {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
            end
            d = obj.distanceMoved(indices);
            if isfield(obj.proc,'rules')
                if strcmp(obj.proc.rules,'diskOn')
                    d(obj.onTarget('diskturn1',indices)) = NaN;
                    d(obj.onTarget('diskturn2',indices)) = NaN;
                    return
                end
            end
            disp('disk turn information unavailable!')
        end

        function [v,c] = getAllVisitCounts(obj,targets)
            % [v,c] = getAllVisitCounts(obj,targets) returns all visit 
            % counts for the given targets rois as struct and cell array.
            % The default target rois are updated for the DoubleCarousel.
            arguments
                obj (1,:)
                targets (1,:) cell = {'diskturn1','diskstop1','diskturn2','diskstop2','border','ring1','ring2','food','empty'}
            end
            [v,c] = getAllVisitCounts@CaroDataHandleList(obj,targets);
        end

    end

end
