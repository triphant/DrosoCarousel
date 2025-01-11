classdef CaroDataShuffler < GenericDataHandleList & CaroInterface
%CARODATASHUFFLER Shuffles the data within a list of DrosoCarousel
%   experiments
%
%   Author: Tilman Triphan, tilman.triphan@uni-leipzig.de
%   License: GNU General Public License v3.0
%   Copyright (C) 2024-2025  Tilman Triphan
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
        cdhl
    end

    methods

        function obj = CaroDataShuffler(cdhl)
            arguments
                cdhl (1,:) CaroDataHandleList
            end
            if numel(unique(cdhl.getLength)) ~= 1
                error("All lists in 'cdhl' have to be of identical length!")
            end
            obj.cdhl = cdhl;
            obj.conf.roilist = obj.cdhl(1).conf.roilist;
            obj.conf.keepTime = true;
            obj.conf.shuffle = [];
            obj.conf.orig = false;
        end

        function cdhl = getCaroDataHandleList(obj)
            % cdhl = getCaroDataHandleList(obj) returns the CaroDataHandleList
            % that is to be shuffled
            cdhl = obj.cdhl;
        end

        function n = getLength(obj)
            % n = getLength(obj) returns the summed length
            n = sum(obj.cdhl.getLength);
        end

        function n = listLength(obj)
            % n = listLength(obj) returns the length of the individual
            % experiments
            n = obj.cdhl(1).getLength;
        end

        function splitlist = splitLists(obj)
            % splitlist = splitLists(obj) splits experiments by videos and
            % creates one list containing all videos
            c = arrayfun(@(x)(x.split), obj.cdhl, "UniformOutput", false);
            splitlist = [c{:}];
        end

        function idx = shuffle(obj,params)
            % idx = shuffle(obj,params) shuffles the videos (i.e. blocks) 
            % between experiments. With 'keepTime' the temporal position of
            % the videos is kept only the experiment changes. You can save
            % and load the shuffled orders
            arguments
                obj (1,:)
                params.keepTime (1,1) logical = obj.conf.keepTime;
                params.saveAs (1,1) string = "";
                params.loadAs (1,1) string = "";
            end
            if exist(params.loadAs,"file") > 0
                idx = readmatrix(params.loadAs);
                obj.conf.shuffle = idx;
                return
            end
            idx = nan(numel(obj.cdhl),obj.listLength);
            if params.keepTime
                for i=1:obj.listLength
                     idx(:,i) = (randperm(numel(obj.cdhl))-1)*obj.listLength+i;
                end
            else
                idx = reshape(randperm(obj.getLength),numel(obj.cdhl),obj.listLength);
            end
            obj.conf.shuffle = idx;
            if (params.saveAs.strlength) > 0
                writematrix(idx,params.saveAs);
            end
        end

        function f = fractionOnTarget(obj,target,reshuffle)
            % f = fractionOnTarget(obj,target,reshuffle) calculates
            % fraction of target with the option to shuffle videos (blocks)
            % around
            arguments
                obj (1,:)
                target {CaroInterface.mustBeRoiName(target,obj)}
                reshuffle (1,1) logical = false
            end
            splitlist = obj.splitLists;
            f_split = splitlist.fractionOnTarget(target);
            if isempty(obj.conf.shuffle) || reshuffle
                idx = obj.shuffle;
            else
                idx = obj.conf.shuffle;
            end
            if obj.conf.orig
                idx = reshape(1:numel(splitlist),obj.listLength,[])';
            end
            f = mean(f_split(idx),2);
        end

        function f = isActive(obj,reshuffle)
            % f = isActive(obj,reshuffle) calculates activity of the flies
            % with the option to shuffle videos (blocks) around
            arguments
                obj (1,:)
                reshuffle (1,1) logical = false
            end
            splitlist = obj.splitLists;
            d = arrayfun(@(x) mean(x.isActive,"omitnan"), splitlist);
            if isempty(obj.conf.shuffle) || reshuffle
                idx = obj.shuffle;
            else
                idx = obj.conf.shuffle;
            end
            f = mean(d(idx),2);
        end

    end
end