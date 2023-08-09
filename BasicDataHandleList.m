classdef BasicDataHandleList < matlab.mixin.Copyable
    %BASICDATAHANDLELIST  Base class to handle experiments.
    %   It provides basic functionality for finding and selecting data
    %
    %   Author: Tilman Triphan, tilman.triphan@uni-leipzig.de
    %   License: GNU General Public License v3.0
    %   Copyright (C) 2016-2023  Tilman Triphan
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

%#ok<*JAPIMATHWORKS>

    properties (SetAccess = protected, GetAccess = protected)
        data % struct array containing data
        info % struct containing meta data and other information
    end

    methods

        function obj = BasicDataHandleList(data)
            if nargin > 0
                obj.data = data;
            end
        end

        function d = getData(obj,indices)
            % d = getData(obj,indices) returns the data struct
            if nargin == 1
                indices = ':';
            end
            d = obj.data(indices);
        end

        function v = getInfo(obj,infofield)
            % v = getInfo(obj,infofield) returns information about the list
            if nargin == 1
                v = [obj.info];
            else
                obj.checkForError('invalidInfoField',{infofield});
                v = arrayfun(@(o) o.info.(infofield),obj(:),'Uniform',0);
            end
        end

        function saveToFile(obj,filename)
            % saveToFile(obj,filename) saves the object
            disp('Saving list, please wait...')
            dhl = obj.copy;
            save(filename,'dhl','-v7.3');
        end

        function n = getLength(obj)
            % n = getLength(obj) returns the size of the wrapped data
            n = cell2mat(arrayfun(@(o) numel(o.data),obj(:),'Uniform',0));
        end

        function obj = discardEmptyLists(obj)
            % obj = discardEmptyLists(obj) removes empty lists, i.e. lists 
            % with no data
            obj(obj.getLength == 0) = [];
        end

        function v = addFieldAndValues(obj,fieldname,values)
            % v = addFieldAndValues(obj,fieldname,values) adds a new field 
            % (and values for it) to the data struct.
            obj.checkForError('inputMismatch',{numel(values)});
            if ~iscell(values)
                values = num2cell(values);
            end
            [obj.data(:).(fieldname)] = deal(values{:});
            v = obj.get(fieldname);
        end

        function removeField(obj,fieldname)
            % removeField(obj,fieldname) removes a field from the data 
            % struct
            for o=1:numel(obj)
                obj(o).checkForError('invalidFieldname',{fieldname});
                obj(o).data = rmfield(obj(o).data,fieldname);
            end
        end

        function v = get(obj,fieldname,indices)
            % v = get(obj,fieldname,indices) is used to get the underlying
            % data. Without parameters, it returns a list of
            % matlab:hyperlinks that can be clicked to get to the
            % respective fields. The parameter fieldname is the name of the
            % field of the data struct or the name of the method that
            % will return an array of results. The optional second
            % parameter specifies indices if only a subset of the data
            % should be returned.
            if nargin == 1
                fns = fieldnames(obj.data);
                hist = com.mathworks.mlservices.MLCommandHistoryServices.getSessionHistory;
                lc = char(hist(end));
                for i=1:numel(fns)
                    cmd = strcat('<a href="matlab: ',lc(1:end-4),...
                        '.get(''',fns(i),''')">',fns(i),'</a>');
                    disp(cell2mat(cmd));
                end
                return
            elseif nargin == 2
                % numeric parameter -> return data for index
                if isnumeric(fieldname) || islogical(fieldname)
                    v = obj.getData(fieldname);
                    return
                else
                    indices = ':';
                end
            end
            if ismethod(obj,fieldname)
                v = obj.(fieldname)(indices);
            elseif isprop(obj,fieldname)
                v = obj.(fieldname);
            else
                fn =  strsplit(fieldname,'.');
                if length(fn) == 1
                    v = {obj.data(indices).(fieldname)}';
                else
                    v = obj.(cell2mat(fn(1))).(cell2mat(fn(2)));
                    v = {v(indices)'};
                end
            end
            if ~iscellstr(v) && iscell(v)
                if isobject(v{1})
                    v = [v{:}]';
                else
                    v = cell2mat(v);
                end
            end
        end

        function dhl = selectListsByInfoQuery(obj,varargin)
            % dhl = selectListsByInfoQuery(obj,varargin) creates a new 
            % array of dhls selected from an array of source lists based on 
            % a query
            dhl = obj(obj.query(varargin{:},'info')).copy;
        end

        function l = query(obj,queryparams,level)
            % l = query(obj,queryparams,level) queries the data struct or 
            % the info struct
            if nargin == 2
                level = 'data';
            end
            n = round(numel(queryparams)/3);
            params = reshape(queryparams,3,n);
            keys = params(1,:);
            comp = params(2,:);
            vals = params(3,:);
            switch level
                case 'info'
                    idxmat = false(numel(obj),n);
                case 'data'
                    idxmat = false(obj.getLength,n);
            end
            for i=1:numel(comp)
                switch level
                    case 'info'
                        res = obj.getInfo(keys{i});
                    case 'data'
                        res = obj.get(keys{i});
                end
                if iscell(res)
                    if isnumeric(res{1})
                        res = cell2mat(res);
                    end
                end
                switch comp{i}
                    case {'=','==','eq','is'}
                        if isnumeric(res)
                            indices = (res == str2double(vals{i}));
                        elseif islogical(res{1})
                            if str2double(vals{i}) == 1
                                indices = cell2mat(res);
                            else
                                indices = ~cell2mat(res);
                            end
                        else
                            indices = strcmp(res,vals{i});
                        end
                    case {'~','~=','!=','ne','not'}
                        if isnumeric(res)
                            indices = (res ~= str2double(vals{i}));
                        elseif islogical(res{1})
                            if str2double(vals{i}) == 1
                                indices = ~cell2mat(res);
                            else
                                indices = cell2mat(res);
                            end
                        else
                            indices = ~strcmp(res,vals{i});
                        end
                    case {'>','gt','greater'}
                        if isnumeric(res)
                            indices = (res > str2double(vals{i}));
                        else
                            warning([comp{i} ' ignored for non-numeric parameter ' keys{i}])
                        end
                    case {'>=','ge'}
                        if isnumeric(res)
                            indices = (res >= str2double(vals{i}));
                        else
                            warning([comp{i} ' ignored for non-numeric parameter ' keys{i}])
                        end
                    case {'<','lt','less'}
                        if isnumeric(res)
                            indices = (res < str2double(vals{i}));
                        else
                            warning([comp{i} ' ignored for non-numeric parameter ' keys{i}])
                        end
                    case {'<=','le'}
                        if isnumeric(res)
                            indices = (res <= str2double(vals{i}));
                        else
                            warning([comp{i} ' ignored for non-numeric parameter ' keys{i}])
                        end
                    case {'><','between'}
                        if isnumeric(res)
                            v = str2num(vals{i}); %#ok<ST2NM>
                            if numel(v) == 2
                                indices = (res > v(1)) & (res < v(2));
                            else
                                warning('Two values needed!')
                            end
                        else
                            warning([comp{i} ' ignored for non-numeric parameter ' keys{i}])
                        end
                    case {'contains','like'}
                        indices = ~cellfun(@isempty,strfind(res,vals{i}));
                    case {'match','regexp'}
                        indices = ~cellfun(@isempty,regexp(res,vals{i}));
                    case {'in','ismember'}
                        indices = ismember(res,vals{i});
                    otherwise
                        disp([comp{i} ' not recognized'])
                end
                idxmat(:,i) = indices;
            end
            l = logical(min(idxmat,[],2));
        end

        function indices = findDataByFieldValue(obj,fieldname,fieldvalue)
            % indices = findDataByFieldValue(obj,fieldname,fieldvalue)
            indices = obj.query({fieldname,'==',fieldvalue});
        end

        function removeDataByIndex(obj,indices)
            % removeDataByIndex(obj,indices) deletes selected data from the
            % underlying struct
            old_length = obj.getLength;
            obj.data(indices) = [];
            sprintf('old length: %d, new length: %d',old_length, obj.getLength)
        end

        function removeDataByFieldValue(obj,fieldname,fieldvalue)
            % removeDataByFieldValue(obj,fieldname,fieldvalue) deletes data
            % from underlying struct for data where specified field has a 
            % given value
            indices = obj.findDataByFieldValue(fieldname,fieldvalue);
            obj.removeDataByIndex(indices);
        end

        function setLength(obj,length)
            % setLength(obj,length) shortens the list to length
            if nargin == 1
                return
            end
            for o=1:numel(obj)
                obj(o).data = obj(o).getData(1:min(length,obj(o).getLength));
            end
        end

        function sortDataByFieldValue(obj,fieldname,mode)
            % sortDataByFieldValue(obj,fieldname,mode) sorts the underlying
            % data by values for a field
            if nargin == 2
                mode = 'ascend';
            end
            for o=1:numel(obj)
                d = obj(o).get(fieldname);
                if iscell(d)
                    [~, order] = sort(d);
                else
                    [~, order] = sort(d,mode);
                end
                obj(o).data = obj(o).data(order);
            end
        end

        function dhl = createSubsetByIndex(obj,indices)
            % dhl = createSubsetByIndex(obj,indices) creates a subset list
            % from given indices
            dhl = obj.copy;
            dhl.data = dhl.data(indices);
        end

        function dhl = subs(obj,indices)
            % dhl = subs(obj,indices) is a shortcut for createSubsetByIndex
            dhl = obj.createSubsetByIndex(indices);
        end

        function dhls = split(obj)
            % dhls = split(obj) splits the list into sublists for each data
            % item
            dhls = feval([class(obj),'.empty'],obj.getLength,0);
            for i = 1:obj.getLength
                dhls(i) = obj.createSubsetByIndex(i);
            end
        end

        function dhl = createSubsetByFieldValue(obj,fieldname,fieldvalue)
            % dhl = createSubsetByFieldValue(obj,fieldname,fieldvalue) only
            % keeps data for a given value in a field
            indices = obj.findDataByFieldValue(fieldname,fieldvalue);
            dhl = obj.createSubsetByIndex(indices);
        end

        function dhl = createSubsetByQuery(obj,queryparams,level)
            % dhl = createSubsetByQuery(obj,queryparams,level) creates a 
            % subset list for data where a certain condition is true
            if nargin == 2
                level = 'data';
            end
            dhl = feval([class(obj),'.empty'],numel(obj),0);
            for o=1:numel(obj)
                indices = obj(o).query(queryparams,level);
                dhl(o) = obj(o).createSubsetByIndex(indices);
            end
        end

        function addInfo(obj,fieldname,values)
            % addInfo(obj,fieldname,values) adds or updates info structure
            if numel(obj) == 1
                obj.info.(fieldname) = values;
            else
                for o=1:numel(obj)
                    obj(o).info.(fieldname) = values{o};
                end
            end
        end

    end

    methods (Access = protected)

        function this = findThis(obj) %#ok<MANU>
            % this = findThis(obj) returns the name of the object variable
            hist = com.mathworks.mlservices.MLCommandHistoryServices.getSessionHistory;
            lastcom = char(hist(end));
            this = lastcom(1:strfind(lastcom,'.')-1);
        end

        function checkForError(obj,errorname,params)
            % checkForError(obj,errorname,params)
            switch errorname
                case 'fileNotFound'
                    if ~exist(params{1},'file')
                        error('BasicDataHandleList:fileNotFound',...
                            'Error. \nFile ''%s'' not found!',params{1})
                    end
                case 'invalidObjectLength'
                    if numel(obj) ~= params{1}
                        error('BasicDataHandleList:invalidObjectLength',...
                            'Error. \nThis method can only be applied on 1x%d lists.',params{1})
                    end
                case 'differentObjectLength'
                    if numel(obj) ~= numel(params{1})
                        error('BasicDataHandleList:differentObjectLength',...
                            'Error. \nThe two objects differ in length.')
                    end
                case 'invalidFieldname'
                    if ~isfield(obj.data,params{1})
                        error('BasicDataHandleList:invalidFieldname',...
                            'Error. \nA field named ''%s'' doesn''t exist!',params{1})
                    end
                case 'duplicateFieldname'
                    if isfield(obj.data,params{1})
                        error('BasicDataHandleList:duplicateFieldname',...
                            'Error. \nA field named ''%s'' already exists!',params{1})
                    end
                case 'invalidInfoField'
                    if ~isfield(obj(1).info,params{1})
                        error('BasicDataHandleList:invalidInfoField',...
                            'Error. \nA infofield named ''%s'' doesn''t exist!',params{1})
                    end
                case 'differentListTypes'
                    if params{1} ~= 1
                        error('BasicDataHandleList:differentListTypes',...
                            'Error. \nMore than one type of list found!')
                    end
                case 'inputMismatch'
                    if ~ismember(params{1}, [1 obj.getLength])
                        error('BasicDataHandleList:inputMismatch',...
                            'Error. \nNumber of values does not match!')
                    end
            end
        end

    end

end
