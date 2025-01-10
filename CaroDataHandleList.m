classdef CaroDataHandleList < BasicDataHandleList & CaroInterface
%CARODATAHANDLELIST Class to handle DrosoCarousel experiments
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

    properties (Dependent, Hidden)
        ActiveRois
        StartTime
        Offset
        MaxDataPoints
    end

    methods

        function obj = CaroDataHandleList(data)
            if nargin > 0
                obj.data = data;
            end
        end

        function v = get.StartTime(obj)
            % v = get.StartTime(obj) returns the start time of the
            % experiment
            c = obj.getInfo('startTime');
            v = c{:};
        end

        function v = get.Offset(obj)
            % v = get.Offset(obj) returns the offset (in frames) for the
            % experiment
            v = seconds(timeofday(obj.StartTime))*obj.ctrl.framesPerSec;
        end

        function v = get.ActiveRois(obj)
            % v = get.ActiveRois(obj) returns the active ROIs
            ac = obj.conf.roiprio;
            v = obj.conf.roilist(ac);
        end

        function d = get.MaxDataPoints(obj)
            % d = get.MaxDataPoints(obj) returns the number of expected
            % frames per 24h window
            d = obj.ctrl.maxFrame*(60/obj.ctrl.videoLength)*24;
        end

        function d = createEthogram(obj,indices,params)
            % d = createEthogram(obj,indices,params) creates a pseudo-
            % ethogram and shows for each frame the ROI the fly is in.
            % Use 'shuffleEvents=true' to shuffle the position of events
            arguments
                obj (1,1)
                indices {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
                params.shuffleEvents (1,1) logical = false;
            end
            rois = obj.ActiveRois;
            d = zeros(numel(obj.get('x',indices)),1);
            for f = 1:numel(rois)
                d(obj.onTarget(rois{f},indices)) = obj.conf.roiprio(f);
            end
            if params.shuffleEvents
                evts = find(d - circshift(d,1) ~= 0);
                n = diff([evts; numel(d)+1]);
                idx = randperm(numel(evts));
                s = [1; cumsum(n(idx))];
                v = nan(numel(d),1);
                for i=1:numel(idx)
                    v(s(i):s(i)+n(idx(i))-1) = d(evts(idx(i)):evts(idx(i))+n(idx(i))-1);
                end
                d = v;
            end
        end

        function d = createEthogramPerDay(obj,params)
            % d = createEthogramPerDay(obj,params) creates a pseudo-
            % ethogram for the active ROIs, aligned to 24h windows.
            % Use 'shuffleEvents=true' to shuffle the position of events
            arguments
                obj (1,1)
                params.shuffleEvents (1,1) logical = false;
            end
            maxDays = obj.maxDays;
            offs = obj.Offset;
            vals = nan(maxDays*obj.MaxDataPoints,1);
            etho = obj.createEthogram('shuffleEvents',params.shuffleEvents);
            vals(offs+1:offs+numel(etho)) = etho;
            %label points outside of the experiment as 'nodata'
            vals(isnan(vals)) = find(strcmp(obj.getRoiList,'nodata'));
            d = reshape(vals,[],maxDays);
        end

        function l = findTransitions(obj,target,indices)
            % l = findTransitions(obj,target,indices) finds transitions
            % between target and other ROIs.
            arguments
                obj (1,1)
                target (1,:) char {CaroInterface.mustBeRoiName(target,obj)}
                indices {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
            end
            ot = obj.onTarget(target,indices);
            l = logical([0; diff(ot)]);
        end

        function [lon,loff,ton,toff,tdur,tint] = getVisitMetrics(obj,target,indices,params)
            % [lon,loff,ton,toff,tdur,tint] = getVisitMetrics(obj,target,indices)
            % returns parameters of ROI transitions.
            arguments
                obj (1,:)
                target (1,:) char {CaroInterface.mustBeRoiName(target,obj)}
                indices {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
                params.minVisitLength (1,1) double = obj.conf.minVisitLength;
                params.roiGapMax (1,1) double = obj.conf.roiGapMax;
            end
            ft = obj.findTransitions(target,indices);
            ot = obj.onTarget(target,indices);
            lon = ft & ot;
            loff = ft & ~ot;
            ton = find(lon);
            toff = find(loff);

            if numel(toff) == 0
                tdur = [];
                tint = NaN;
                return
            end

            if toff(1) == 1
                toff(1) = [];
            end

            if ot(end) == 1
                ton(end) = [];
            end

            if numel(ton) < numel(toff)
                ton = [1; ton];
            end

            tint = circshift(ton,-1)-toff;
            idx = find(tint <= params.roiGapMax);
            
            if numel(tint) == 0
                tdur = [];
                tint = NaN;
                return
            end

            idx(end) = [];
            toff(idx) = [];
            ton(idx+1) = [];

            tdur = toff-ton;
            idx = tdur < params.minVisitLength;
            ton(idx) = [];
            toff(idx) = [];

            tdur = toff-ton;
            tint = circshift(ton,-1)-toff;
            if ~isempty(tint)
                tint(end) = [];
            end
        end

        function stats = findSpecialTransitions(obj,target,params)
            % stats = findSpecialTransitions(obj,target,params) finds 
            % transition events between ROIs with a high rate of occupancy
            % on a second ROI in a given time window.
            arguments
                obj (1,:)
                target (1,:) char {CaroInterface.mustBeRoiName(target,obj)}
                params.indices  {CaroInterface.mustBeValidIndex(params.indices,obj)} = ':'
                params.drawingwindow (1,1) {mustBeInteger} = 500
                params.introsize (1,1) {mustBeInteger} = 50
                params.maxdrawevents (1,1) {mustBeInteger} = 5
                params.mincutoff (1,1) {mustBeInteger} = obj(1).conf.minVisitLength
                params.minfraction (1,1) {mustBeInRange(params.minfraction,0,1)} = 0.25
                params.nofigure (1,1) logical = false
                params.sortby (1,:) char {CaroInterface.mustBeRoiName(params.sortby,obj)}
                params.windowsize (1,1) {mustBeInteger} = 500
            end
            for o=1:numel(obj)
                al = obj(o).conf.roilist(obj(o).conf.roicode);
                sbCode = ismember(al,params.sortby);
                for i=1:numel(al)
                    a = al(i);
                    lon.(cell2mat(a)) = obj(o).onTarget(cell2mat(a));
                end
                [~,~,~,toff,tdur,~] = obj(o).getVisitMetrics(target,params.indices);
                remdat = tdur < params.mincutoff;
                toff(remdat) = [];
                tdur(remdat) = [];
                fns = fieldnames(lon);
                f = [];
                dl = obj(o).getDataLength;
                for i = 1:numel(toff)
                    if numel(toff) > 0
                        for fn = 1:numel(fns)
                            l2 = lon.(fns{fn});
                            ep = min(toff(i)+params.windowsize,dl);
                            f(i).(fns{fn}) = sum(l2(toff(i):ep))/min(params.windowsize,ep-toff(i)); %#ok<AGROW> 
                        end
                    end
                    for fn = 1:numel(fns)
                        f.(fns{fn}) = NaN;
                    end
                end
                frc = nan(numel(tdur),numel(fns));
                for fn = 1:numel(fns)
                    frc(:,fn) = [f.(fns{fn})];
                end
                [~,idx] = sort(frc(:,sbCode),'desc');
                frc = frc(idx,:);
                if ~params.nofigure
                    h = figure;
                    hold on
                    bar(frc,'stacked');
                    ac = obj(o).conf.roicode;
                    for i=1:numel(h.CurrentAxes.Children)
                        c = numel(h.CurrentAxes.Children)-i+1;
                        h.CurrentAxes.Children(i).FaceColor = obj(o).conf.cmapdefault(ac(c),:);
                    end
                    xlabel('Events')
                    ylabel('Fractions')
                    ylim([0 1])
                    legend(fns)
                end
                eidx = frc(:,sbCode) >= params.minfraction;
                stats(o).count = sum(eidx);
                stats(o).frac = frc(eidx,sbCode);
                stats(o).relfrac = sum(eidx)/numel(eidx);
                if params.maxdrawevents > 0 && ~params.nofigure
                    temp_im = obj(o).conf.indexMode;
                    obj(o).conf.indexMode = 'frame';
                    for z=1:min(params.maxdrawevents,numel(idx))
                        fs = max(toff(idx(z))-params.drawingwindow,1);
                        fe = min(toff(idx(z))+params.drawingwindow,obj(o).getDataLength);
                        stats(o).h(z+1) = obj(o).plotWalkingTrace(fs:fe);
                    end
                    obj(o).conf.indexMode = temp_im;
                end
                stats(o).fidx = toff(idx(eidx));
            end
        end

        function h = plotVisitDurationHistogram(obj,target,indices,params)
            % h = plotVisitDurationHistogram(obj,target,indices,params) 
            % plots a histogram of the duration of visits to a ROI.
            arguments
                obj (1,:)
                target (1,:) char {CaroInterface.mustBeRoiName(target,obj)}
                indices {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
                params.edges (1,:) char {mustBeMember(params.edges,{'log10','log2'})} = 'log10'
            end
            h = figure;
            h.Name = ['Visit durations on ',target];
            if strcmp(params.edges,'log10')
                edges = logspace(0,6,7);
            else
                edges = arrayfun(@(x) 2^x,1:20);
            end
            d = nan(numel(obj),numel(edges)-1);
            for o=1:numel(obj)
                [~,~,~,~,tdur] = obj(o).getVisitMetrics(target,indices);
                d(o,:) = histcounts(tdur,edges);
            end
            if numel(obj) == 1
                bar(d)
            else
                bar(mean(d))
            end
            xlabel('Duration of visit [frames]')
            xticks(0.5:1:numel(edges))
            xticklabels(edges)
            ylabel('Number of visits (on)')
        end

        function cdhls = shortenDataToTimeWindow(obj,params)
            % cdhls = shortenDataToTimeWindow(obj,params) returns a new 
            % list or lists with the data cropped to a specific time window
            arguments
                obj (1,:)
                params.startTime (1,1) double {mustBeInteger,mustBeInRange(params.startTime,0,24)} = 6
                params.duration (1,1) double {mustBePositive,mustBeInteger} = 16
                params.offset (1,1) double {mustBeInteger} = 0
                params.endHour (1,1) double {mustBeInteger,mustBeInRange(params.endHour,-1,24)} = -1
                params.endDay (1,1) double {mustBeInteger} = 0
            end
            cdhls = feval([class(obj),'.empty'],numel(obj),0);
            for o=1:numel(obj)
                if params.endDay == 0
                    si = find(obj(o).get('dateTime').Hour == params.startTime,1);
                    idx = si+params.offset*2:si+params.duration*2-1+params.offset*2;
                    if obj(o).getLength < si+params.duration*2-1+params.offset*2
                        continue
                    end
                    cdhls(o) = obj(o).createSubsetByIndex(idx);
                    cdhls(o).info.startTime = cdhls(o).get('dateTime',1);
                else
                    dt = obj(o).get('dateTime');
                    dnum = floor(datenum(dt));
                    dn = dnum - dnum(1);
                    idx = find(dn == params.endDay, 1, 'last');
                    if params.endHour ~= -1
                        id = dn == params.endDay;
                        ih = dt.Hour == params.endHour;
                        idx = find(all([id ih],2),1,'first')-1;
                    end
                    if ~isempty(idx)
                        cdhls(o) = obj(o).createSubsetByIndex(1:idx);
                    end
                end
            end
            cdhls = cdhls.discardEmptyLists;
            fprintf('%d of %d fit into window\n',numel(cdhls),numel(obj))
        end

        function h = plotWalkingTrace(obj,indices)
            % h = plotWalkingTrace(obj,indices) plots a 3D walking trace 
            % with time as the 3rd dimension and a temporal color code
            arguments
                obj (1,1)
                indices  {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
            end
            h = gobjects(numel(obj),1);
            for o = 1:numel(obj)
                h(o) = figure;
                h(o).Name = ['WalkingTrace for ',obj(o).info.name];
                x = obj(o).get('x',indices);
                y = obj(o).get('y',indices);
                c = 1:numel(x);
                surface([x,x], [y,y], [c(:),c(:)], 'EdgeColor','flat', 'FaceColor','none');
                colormap(jet(numel(x)));
                set(gca,'Ydir','reverse');
                zlabel('Time');
                axis([0 obj(o).ctrl.imageResX 0 obj(o).ctrl.imageResY]);
                axis square
            end
        end

        function h = plotPositions(obj,indices,params)
            % h = plotPositions(obj,indices,params) plots x/y positions
            % with time on the z-axis
            arguments
                obj (1,1)
                indices  {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
                params.plot (1,:) char {mustBeMember(params.plot,{'roi','time'})} = 'roi'
                params.size (1,1) double = 3
            end
            disp('Creating plot(s)...')

            h = gobjects(numel(obj),1);
            for o=1:numel(obj)
                h(o) = figure;
                x = obj(o).get('x',indices);
                y = obj(o).get('y',indices);
                switch params.plot
                    case 'roi'
                        h(o).Name = ['Positions (ROI coded) for ',obj(o).info.name];
                        fns = obj(o).conf.roilist(obj(o).conf.roiprio);
                        c = zeros(numel(x),1);
                        for f = 1:numel(fns)
                            c(obj(o).onTarget(fns{f},indices)) = obj(o).conf.roiprio(f);
                        end
                        s = params.size;
                        colormap(obj(o).conf.cmapdefault);
                    case 'time'
                        h(o).Name = ['Positions (time coded) for ',obj(o).info.name];
                        c = 1:numel(x);
                        s = params.size;
                        colormap(jet(numel(x)));
                end
                scatter3(x,y,1:numel(x),s,c,'filled');
                if strcmp(params.plot,'roi')
                    caxis([1 size(obj(o).conf.cmapdefault,1)])
                end
                set(gca,'Ydir','reverse');
                zlabel('Time');
                axis([0 obj(o).ctrl.imageResX 0 obj(o).ctrl.imageResY]);
                axis square
                view(45,30)
            end
        end

        function h = plotOccupancy(obj,indices)
            % h = plotOccupancy(obj,indices) plots a 2D histogram of 
            % occupancy
            arguments
                obj (1,1)
                indices  {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
            end
            h = gobjects(numel(obj),1);
            for o = 1:numel(obj)
                h(o) = figure;
                h(o).Name = ['Occupancy for ',obj(o).info.name];
                x = obj(o).get('x',indices);
                y = obj(o).get('y',indices);
                resX = obj(o).ctrl.imageResX*obj(o).conf.occupancyScale;
                resY = obj(o).ctrl.imageResY*obj(o).conf.occupancyScale;
                image(hist3([[0;x(:);obj(o).ctrl.imageResX], ...
                    [0;y(:);obj(o).ctrl.imageResY]], [resX resY])');
                colormap(hot);
                axis square
            end
        end

        function h = plotMeanOccupancy(obj,indices,params)
            % h = plotMeanOccupancy(obj,indices,params) plots a 2D
            % histogram of mean occupancy for several experiments
            arguments
                obj (1,1)
                indices  {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
                params.occupancyScale (1,1) = obj(1).conf.occupancyScale;
            end
            h = figure;
            h.Name = 'Mean occupancy';
            imgResX = unique(arrayfun(@(o) o.ctrl.imageResX,obj));
            imgResY = unique(arrayfun(@(o) o.ctrl.imageResY,obj));
            if numel(imgResX) + numel(imgResY) ~= 2
                error('CaroDataHandleList:inconsistentResolution',...
                    'image resolution for x and y should be consistent')
            end
            resX = imgResX*params.occupancyScale;
            resY = imgResY*params.occupancyScale;
            img = zeros(resX,resY,numel(obj));
            for o=1:numel(obj)
                x = obj(o).get('x',indices) - obj(o).ctrl.rois.arena.x + obj(o).ctrl.imageResX/2;
                y = obj(o).get('y',indices) - obj(o).ctrl.rois.arena.y + obj(o).ctrl.imageResY/2;
                img(:,:,o) = hist3([[0;x(:);obj(o).ctrl.imageResX], ...
                    [0;y(:);obj(o).ctrl.imageResY]], [resX resY])';
            end
            image(mean(img,3));
            colormap(hot);
            axis square
        end

        function h = plotEthogram(obj,indices,params)
            % h = plotEthogram(obj,indices,params) plots a pseudo-ethogram
            % i.e. color codes each frame with the respective ROI the fly
            % is in
            arguments
                obj (1,:)
                indices {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
                params.setXAxis (1,1) char {mustBeMember(params.setXAxis,{'f','v','m','h','d'})}= 'f'
                params.shuffleEvents (1,1) logical = false;
            end
            h = gobjects(numel(obj),1);
            for o=1:numel(obj)
                h(o) = figure;
                h(o).Name = ['Ethogram for ',obj(o).info.name];
                cdata = obj(o).createEthogram(indices,'shuffleEvents',params.shuffleEvents)';
                image(cdata);
                colormap(obj(o).conf.cmapdefault);
                caxis([1 size(obj(o).conf.cmapdefault,1)])
                yticklabels('');
                set(gca, 'TickDir', 'out');
                obj(o).addRoiLegend(h(o));
                obj(o).setXAxis(h(o),'Unit',params.setXAxis,'indices',indices);
            end
        end

        function h = plotEthogramPerDay(obj,params)
            % h = plotEthogramPerDay(obj,params) plots a pseudo-ethogram
            % broken down into 24h windows
            arguments
                obj (1,:)
                params.shuffleEvents (1,1) logical = false;
                params.showLegend (1,1) logical = true
                params.legendLocation (1,:) char = 'NorthWest'
            end
            h = gobjects(numel(obj),1);
            for o=1:numel(obj)
                h(o) = figure;
                h(o).Name = ['Ethogram (split by days) for ',obj(o).info.name];
                cdata = obj(o).createEthogramPerDay('shuffleEvents',params.shuffleEvents)';
                ndays = size(cdata,1);
                mdp = obj(o).MaxDataPoints;
                image(cdata);
                colormap(obj(o).conf.cmapdefault);
                yticklabels('');
                set(gca, 'TickDir', 'out');
                axis([0 mdp 0.5 ndays+0.5])
                yticks(0:ndays);
                yticklabels(0:ndays);
                ylabel('days');
                xticks([0 mdp/4 mdp/2 mdp*0.75 mdp*22/24 mdp]);
                xticklabels({'0:00','06:00','12:00','18:00','22:00','24:00'});
                if params.showLegend
                    obj(o).addRoiLegend(h(o),'Location',params.legendLocation);
                end
            end
        end

        function h = plotEthogramPerDays(obj,params)
            % h = plotEthogramPerDays(obj,params) plots a single pseudo-
            % ethogram for several experiments aligned by time
            arguments
                obj (1,:)
                params.titleStr char = ['Ethograms for ',obj.findThis]
                params.makeAvg logical = false
                params.includeRoi cell = obj(1).getActiveRois
                params.excludeRoi cell = {}
                params.moveMeanWin (1,1) double = 1500;
                params.shuffleEvents (1,1) logical = false;
                params.xLimAt (1,1) double = 24;
                params.showActivity (1,1) logical = false;
            end
            h = figure;
            if params.makeAvg
                subplot('Position',[0.1 0.3 0.8 0.65])
            end
            h.Name = params.titleStr;
            maxDays = obj.maxDays;
            mdp = obj(1).MaxDataPoints;
            cdata_all = nan(numel(obj),maxDays*mdp);
            for o=1:numel(obj)
                cdata = reshape(obj(o).createEthogramPerDay('shuffleEvents',params.shuffleEvents),[],1)';
                cdata_all(o,1:numel(cdata)) = cdata;
            end
            cdata_all(isnan(cdata_all)) = find(strcmp(obj(o).getRoiList,'nodata'));
            image(cdata_all);
            colormap(obj(o).conf.cmapdefault);
            title(params.titleStr,'Interpreter','none');
            set(gca, 'TickDir', 'out');
            offset = mdp * 0.25;
            axis([offset maxDays*mdp 0.5 numel(obj)+0.5])
            xt = offset:mdp/2:maxDays*mdp;
            xt(2:2:end) = xt(2:2:end) + mdp / 6;
            xticks(xt);
            xticklabels(repmat({'06:00';'22:00'},maxDays,1));
            yticks(1:numel(obj)+1);
            yticklabels(obj.getInfo('name'));
            set(gca,'TickLabelInterpreter','none');
            obj(o).addRoiLegend(h);
            xend = size(cdata_all,2)-obj(1).ctrl.maxFrame*(1440/obj(1).ctrl.videoLength-params.xLimAt*2);
            xlim([offset xend]);
            if params.makeAvg
                subplot('Position',[0.1 0.05 0.8 0.2])
                hold on
                idx = ~ismember(params.includeRoi,params.excludeRoi);
                al = obj(o).conf.roilist;
                rois = find(ismember(al,params.includeRoi(idx)));
                for i=1:numel(rois)
                    ac = rois(i);
                    roi = obj(o).ctrl.rois.(al{ac});
                    v = sum(cdata_all==ac)/numel(obj);
                    if params.moveMeanWin > 0
                        v = movmean(v,params.moveMeanWin);
                    end
                    plot(v,'Color',roi.color,'LineWidth',1);
                end
                if params.showActivity
                    v = mean(obj.activitySync,2,'omitnan');
                    if params.moveMeanWin > 0
                        v = movmean(v,params.moveMeanWin);
                    end
                    plot(v,'Color',[0 1 1],'LineWidth',1);
                end
                xlim([offset xend]);
                ylim([0 1]);
                xticks(xt);
                xticklabels(repmat({'06:00';'22:00'},maxDays,1));
            end
        end

        function [h,t] = plotFractions(obj,indices,params)
            % [h,t] = plotFractions(obj,indices,params) plots the fractions
            % of time the fly spends in different target ROIs
            arguments
                obj (1,:)
                indices {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
                params.saveTable (1,1) logical = false
            end
            h = figure;
            h.Name = ['Fractions of time in ROI for ',obj.findThis];
            ac = obj(1).conf.roicode;
            fns = obj(1).conf.roilist(ac);
            frc = nan(numel(obj),numel(fns));
            for f = 1:numel(fns)
                frc(:,f) = obj.fractionOnTarget(fns{f},indices);
            end
            if size(frc,1) == 1
                frc = [frc;frc];
            end
            boxplot(frc,'colors',obj(1).conf.cmapdefault(ac,:),...
                'OutlierSize',4,'Symbol','o','labels',fns)
            if params.saveTable
                t = array2table(frc,'VariableNames',fns);
                try
                    writetable(t,[pwd '\' findThis(obj) '_fractions.xlsx']);
                catch
                    disp('failed to write table!')
                end
            else
                t = [];
            end
        end

        function h = plotFractionsStacked(obj,indices)
            % h = plotFractionsStacked(obj,indices) plots the fractions of
            % time the fly spends in different ROIs as stacked bars for 
            % each experiment
            arguments
                obj (1,:)
                indices {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
            end
            h = figure;
            h.Name = ['Fractions of time in ROI for ',obj.findThis];
            ac = obj(1).conf.roicode;
            fns = obj(1).conf.roilist(ac);
            frc = zeros(numel(obj),numel(fns));
            for f = 1:numel(fns)
                frc(:,f) = obj.fractionOnTarget(fns{f},indices);
            end
            if size(frc,1) == 1
                frc = [frc;frc];
            end
            colororder(h,obj(1).conf.cmapdefault(ac,:))
            bar(frc,'stacked');
            axis([0.5 numel(obj)+0.5 0 1])
            xticklabels(obj.getInfo('name'));
            set(gca,'TickLabelInterpreter','none');
            legend(fns)
        end

        function h = plotOnTarget(obj,target,indices,params)
            % h = plotOnTarget(obj,target,indices,params) shows the 
            % distribution of stays on ROI 'target' over the course of the
            % experiment and a percentage value. If 'allRois' is set to 
            % 'true', all ROIs are shown for a single experiment instead
            arguments
                obj (1,:)
                target (1,:) char {CaroInterface.mustBeRoiName(target,obj)}
                indices {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
                params.allRois (1,1) logical = false
            end
            disp('new')
            h = figure;
            if params.allRois && (numel(obj)==1)
                al = obj.conf.roilist(obj.conf.roicode);
                n = numel(al);
            else
                n = numel(obj);
            end
            for o=1:n
                if params.allRois
                    idx = 1;
                    target = al{o};
                else
                    idx = o;
                end
                ax = subplot('position',[0.05, 1-o/(n+2)-0.03, 0.75, 1/(n+10)]);
                image(obj(idx).onTarget(target,indices)');
                yticklabels('');
                colormap(ax,[1 1 1; obj(idx).ctrl.rois.(target).color])
                xlabel(target)
                fn = max(obj(idx).fractionOnTarget(target,indices),10^-10);
                ax = subplot('position',[0.85, 1-o/(n+2)-0.03, 0.2, 1/(n+10)]);
                pie(fn);
                colormap(ax,[1 1 1; obj(idx).ctrl.rois.(target).color])
            end
        end

        function [h,stats] = plotXYPlots(obj,cdhls_array,params)
            % [h,stats] = plotXYPlots(obj,cdhls_array,params) plots X-Y
            % plots for several groups of experiments
            arguments
                obj (1,:)
                cdhls_array (1,:) cell
                params.dataX (1,:) char {CaroInterface.mustBeValidInfoFieldOrMethod(params.dataX,obj)}
                params.dataY (1,:) char {CaroInterface.mustBeValidInfoFieldOrMethod(params.dataY,obj)}
                params.paramX = ''
                params.paramY = ''
                params.legend (1,:) char {CaroInterface.mustBeValidInfoField(params.legend,obj)} = 'genotype'
            end
            disp('Processing, please wait...')
            h = figure;
            hold on
            fns = cell(numel(cdhls_array),1);
            stats = [];
            if ismethod(obj,params.dataX)
                dataXMode = 'method';
            end
            if ismethod(obj,params.dataY)
                dataYMode = 'method';
            end
            if any(strcmp(params.dataX,fieldnames(obj.getInfo)))
                dataXMode = 'info';
            end
            if any(strcmp(params.dataY,fieldnames(obj.getInfo)))
                dataYMode = 'info';
            end
            for a=1:numel(cdhls_array)
                label = cdhls_array{a}.getInfo(params.legend);
                fns{a} = [label{1},' [',num2str(numel(cdhls_array{a})),']'];
                switch dataXMode
                    case 'method'
                        if isempty(params.paramX)
                            xd = cdhls_array{a}.(params.dataX);
                        else
                            xd = cdhls_array{a}.(params.dataX)((params.paramX));
                        end
                    case 'info'
                        xd = cell2mat(cdhls_array{a}.getInfo(params.dataX));
                end
                switch dataYMode
                    case 'method'
                        if isempty(params.paramY)
                            yd = cdhls_array{a}.(params.dataY);
                        else
                            yd = cdhls_array{a}.(params.dataY)((params.paramY));
                        end
                    case 'info'
                        yd = cell2mat(cdhls_array{a}.getInfo(params.dataY));
                end
                s = scatter(xd,yd,30,[bitget(a,1) bitget(a,2) bitget(a,3)],'filled');
                [rho,pval] = corr(xd,yd);
                hasbehavior(s,'legend',false);

                stats(a).group = label;
                stats(a).x = xd;
                stats(a).y = yd;
                stats(a).rho = rho;
                stats(a).pval = pval;
                
                f = polyval(polyfit(xd,yd,1),xd);
                p = plot(xd,f,'color',[bitget(a,1) bitget(a,2) bitget(a,3)]);
                hasbehavior(p,'legend',false);

                xdm = mean(xd,'omitnan');
                ydm = mean(yd,'omitnan');
                xds = std(xd,'omitnan')/sqrt(numel(cdhls_array{a}));
                yds = std(yd,'omitnan')/sqrt(numel(cdhls_array{a}));
                e = errorbar(xdm,ydm,yds,yds,xds,xds,'linestyle','none','color','k');
                scatter(xdm,ydm,100,[bitget(a,1) bitget(a,2) bitget(a,3)],'filled');
                hasbehavior(e,'legend',false);
            end
            xlabel([params.dataX ' ' params.paramX]);
            ylabel([params.dataY ' ' params.paramY]);
            set(gca,'linewidth',1.5)
            set(gca,'fontsize',12)
            legend(fns,'Interpreter','none','Location','best')
        end

        function addFixFile(obj)
            % addFixFile(obj) adds fix-files (i.e. a file containing 
            % manually fixed positions for missed detections) to each
            % experiment (if the file is available)
            arguments
                obj (1,:)
            end
            for o=1:numel(obj)
                if ~isfield(obj(o).info,'fixFile') || isempty(obj(o).info.fixFile)
                    obj(o).info.fixFile = [obj(o).info.dataPath '\' obj(o).info.fullName '.fix'];
                end
                ffn = obj(o).info.fixFile;
                if exist(ffn,'file')
                    obj(o).proc.ri = obj(o).loadFixFile(ffn);
                    obj(o).info.hasFF = true;
                    disp([cell2mat(obj(o).getInfo('name')) ': Fix-file assigned!'])
                else
                    obj(o).info.hasFF = false;
                    disp([cell2mat(obj(o).getInfo('name')) ': No fix-file found!'])
                end
                obj(o).proc.fixApplied = false;
            end
        end

        function applyFixFile(obj,params)
            % applyFixFile(obj,params) applies the fix-file and updates
            % repaired the x and y values
            arguments
                obj (1,:)
                params.version (1,1) double = 1.3
            end
            for o=1:numel(obj)
                if ~isfield(obj(o).proc,'ri')
                    disp([cell2mat(obj(o).getInfo('name')) ': No fix-file assigned!'])
                    continue
                end
                if obj(o).proc.ri.version ~= params.version
                    disp([cell2mat(obj(o).getInfo('name')) ': Invalid fix-file version!'])
                    continue
                end
                for i=1:numel(obj(o).proc.ri.videos)
                    if isempty(obj(o).proc.ri.videos(i).frames)
                        disp([cell2mat(obj(o).getInfo('name')) ': No frames in fix-file, skipping!'])
                        continue
                    end
                    try
                        vn = obj(o).proc.ri.videos(i).fileID;
                        vi = find(obj(o).query({'fileID','eq',vn}));
                        fi = [obj(o).proc.ri.videos(i).frames.index];
                        obj(o).data(vi).x(fi) = [obj(o).proc.ri.videos(i).frames.x];
                        obj(o).data(vi).y(fi) = [obj(o).proc.ri.videos(i).frames.y];
                        obj(o).proc.fixApplied = true;
                    catch
                        disp([cell2mat(obj(o).getInfo('name')) ': Error in fix-file, skipping!'])
                    end
                end
            end
        end

    end

    methods %override

        function v = get(obj,fieldname,indices,indexMode)
            % v = get(obj,fieldname,indices,indexMode) adds an additional
            % indexMode to the get method. You can either adress frames or
            % videos as a whole.
            if nargin < 2
                get@BasicDataHandleList(obj);
                return
            end
            if nargin < 3
                v = get@BasicDataHandleList(obj,fieldname);
                return
            end
            if nargin < 4
                indexMode = obj.conf.indexMode;
            end
            switch indexMode
                case 'video'
                    v = get@BasicDataHandleList(obj,fieldname,indices);
                case 'frame'
                    v = get@BasicDataHandleList(obj,fieldname);
                    v = v(indices);
                otherwise
                    disp("indexMode has to be either 'video' or 'frame'")
            end
        end

    end

    methods (Access = protected)

        function this = findThis(obj,methodname)
            % this = findThis(obj,methodname) returns the name of the
            % object variable
            if nargin < 2
                this = findThis@BasicDataHandleList(obj);
            else
                hist = com.mathworks.mlservices.MLCommandHistoryServices.getSessionHistory; %#ok<JAPIMATHWORKS> 
                lastcom = char(hist(end));
                this = lastcom(1:strfind(lastcom,methodname)-2);
            end
        end

        function d = maxDays(obj)
            % d = maxDays(obj) returns the number of days the experiment
            % lasted
            d = max(arrayfun(@(x) numel(unique(string(x.get('dateTime'),'MM_dd'))),obj));
        end

        function d = activitySync(obj,indices)
            % d = activitySync(obj,indices) lines up activity 
            arguments
                obj (1,:)
                indices {CaroInterface.mustBeValidIndex(indices,obj)} = ':'
            end
            maxDays = obj.maxDays;
            mdp = obj(1).MaxDataPoints;
            d = nan(maxDays*mdp,numel(obj));
            for o=1:numel(obj)
                offs = obj(o).Offset;
                vals = nan(maxDays*mdp,1);
                act = obj(o).isActive(indices);
                vals(offs+1:offs+numel(act)) = act;
                d(:,o) = vals;
            end
        end

        function ri = loadFixFile(~,ffn)
            % ri = loadFixFile(~,ffn) loads a fix-file
            arguments
                ~
                ffn char {mustBeFile}
            end
            ri = jsondecode(fileread(ffn));
            for i=1:numel(ri.videos)
                for j=1:numel(ri.videos(i).frames)
                    % Take care of the 'null' values
                    if isempty(ri.videos(i).frames(j).x)
                        ri.videos(i).frames(j).x = NaN;
                        ri.videos(i).frames(j).y = NaN;
                    end
                end
            end
        end

        function setXAxis(obj,h,params)
            % setXAxis(obj,h,params) updates the x-Axis of the plot with
            % handle h
            arguments
                obj (1,1)
                h (1,1) matlab.ui.Figure
                params.Unit (1,:) char {mustBeMember(params.Unit,{'f','m','h','d','v'})} = 'f'
                params.indices {CaroInterface.mustBeValidIndex(params.indices,obj)} = ':'
            end
            switch params.Unit
                case 'f'
                    h.CurrentAxes.XLabel.String = '[frames]';
                    return
                case 'v'
                    xt = cumsum(obj.fc(params.indices));
                    h.CurrentAxes.XTick = xt;
                    h.CurrentAxes.XTickLabel = arrayfun(@(x)num2str(x),1:numel(xt),'UniformOutput',false);
                    h.CurrentAxes.XLabel.String = '[videos]';
                    return
                case 'm'
                    step = obj.ctrl.framesPerSec*60;
                    h.CurrentAxes.XLabel.String = '[min]';
                case 'h'
                    step = obj.ctrl.framesPerSec*60*60;
                    h.CurrentAxes.XLabel.String = '[hours]';
                case 'd'
                    step = obj.ctrl.framesPerSec*60*60*24;
                    h.CurrentAxes.XLabel.String = '[days]';
            end
            xt = 0:step:h.CurrentAxes.XLim(2);
            h.CurrentAxes.XTick = xt;
            h.CurrentAxes.XTickLabel = arrayfun(@(x)num2str(x),0:numel(xt),'UniformOutput',false);
        end

        function addRoiLegend(obj,h,opts)
            % addRoiLegend(obj,h,opts) adds a legend for the ROI color 
            % codes to the plot with handle h
            arguments
                obj(1,1)
                h (1,1) matlab.ui.Figure = gcf
                opts.Location = 'NorthEast'
            end
            ac = obj.conf.roicode;
            al = obj.conf.roilist;
            nv = nan(numel(al(ac)));
            figure(h);
            colororder(h,obj.conf.cmapdefault(ac,:));
            hold on
            bar(nv,'stacked');
            legend(al(ac),'Location',opts.Location);
        end


    end

end