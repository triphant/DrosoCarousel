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

    end

end
