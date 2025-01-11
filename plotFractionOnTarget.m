function h = plotFractionOnTarget(groups,params)
%PLOTFRACTIONONTARGET Plot fraction on target for groups with the option to
%plot shuffled data
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

    arguments
        groups (1,:) cell
        params.target (1,1) string = "disk"
        params.legend (1,:) string
        params.center (1,1) string {mustBeMember(params.center,["mean","median"])} = "median"
        params.op (1,1) string {mustBeMember(params.op,["divide","subtract"])} = "divide"
        params.scale (1,1) string {mustBeMember(params.scale,["none","log","ln"])} = "none"
        params.abs (1,1) logical = true
    end
    h = figure;
    hold on
    for g=1:numel(groups)
        cdhl = groups{g};
        
        fd = sort(cdhl.fractionOnTarget(params.target),'descend');
        switch params.center
            case "mean"
                fd_c = mean(fd,"omitnan");
            case "median"
                fd_c = median(fd,"omitnan");
        end
        switch params.op
            case "divide"
                fd_op = fd / fd_c;
            case "subtract"
                fd_op = fd - fd_c;
        end
        switch params.scale
            case "none"
                fd_op_s = fd_op;
            case "log"
                fd_op_s = log(fd_op);
            case "ln"
                fd_op_s = log2(fd_op);
        end
        if params.abs
            fd_op_s = abs(fd_op_s);
        end
        fd_op_s(isinf(fd_op_s)) = 15; % (handle infinity values)
        switch class(cdhl)
            case "CaroDataHandleList"
                scatter(-numel(cdhl)/2+1:numel(cdhl)/2,fd_op_s);
            case "CaroDataShuffler"
                n = numel(cdhl.cdhl);
                scatter(-n/2+1:n/2,fd_op_s);
        end
    end
    legend(params.legend)
end