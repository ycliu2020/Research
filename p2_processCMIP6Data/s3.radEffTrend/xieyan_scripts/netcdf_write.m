function [ status ] = netcdf_write(filename,data,structure,if_disp)
% function [ status ] = netcdf_write(filename,data,structure,if_disp)
%
% ARGIN:
%   filename 
%   data
%   structure.name - var/dim/att name string 
%             type - 'dim' / 'att' / 'var'
%             nc_type - 'NC_BYTE' / 'NC_CHAR' / 'NC_SHORT' / 'NC_INT' / 
%                       'NC_FLOAT' / 'NC_DOUBLE'
%             att - structure, e.g. att.mising_value -9999
%             dim - if numeric, prod(dim) must equal length(data(:))
%                   if cellstr, must be existent dimensions
%             var_name - used only when type equals to 'att';
%                        if specified, variable attribute; 
%                        otherwise, global attribute 
%             if_unlim - used only when type equals to 'dim';
%                        if 1, set as NC_UNLIMITED dimension  
%                   
% ARGOUT
%   status - 0: success; -1: error 
%
% note: this program is based on the Matlab netCDF C Interface, which
%   currently doesn't seem to work for the unlimited dimension. (yih.090130)
   
%clear all; close all;
%filename = 'test.nc';
%data = rand(10,10);
%structure = struct('name','rand','type','var','nc_type','NC_DOUBLE','var_name',[]);
% [ status ] = netcdf_write(filename,data,structure)

if ~exist('if_disp','var') || isempty(if_disp)
  if_disp = 0;
end

if exist(filename,'file')
  if if_disp, disp(['Appending to - ',filename]); end
  nc = netcdf.open(filename,'NC_WRITE');
  netcdf.reDef(nc);
else
  if if_disp, disp(['Creating - ',filename]); end
  nc = netcdf.create(filename,'NC_WRITE');
end

if ~exist('structure','var') 
  structure.name = inputname(2);
  structure.type = 'var';
end

if ~isfield(structure,'name')
  structure.name = inputname(2);
end
if ~isfield(structure,'type')
  structure.type = 'var';
end
if ~isfield(structure,'nc_type')
  if isinteger(data)
    structure.nc_type = 'NC_INT';
  elseif isfloat(data)
    structure.nc_type = 'NC_DOUBLE';
  elseif isnumeric(data)
    structure.nc_type = 'NC_DOUBLE';
  elseif isstr(data)
    structure.nc_type = 'NC_CHAR';
  elseif islogical(data)
    structure.nc_type = 'NC_BYTE';
  end
end
if strcmpi(structure.nc_type,'NC_BYTE')
  data = int8(data);
elseif strcmpi(structure.nc_type,'NC_CHAR')
  data = char(data);
elseif strcmpi(structure.nc_type,'NC_SHORT')
  data = int16(data);
elseif strcmpi(structure.nc_type,'NC_INT')
  data = int32(data);
elseif strcmpi(structure.nc_type,'NC_FLOAT')
  data = single(data);
elseif strcmpi(structure.nc_type,'NC_DOUBLE')
  data = double(data);
else
  disp(['ERROR: unrecognizable nc_type!']);
  netcdf.close(nc);
  status = -1;
  return;
end

if strcmpi(structure.type,'variable') ||  strcmpi(structure.type,'var')

  % dimension
  if isfield(structure,'dim') 
    nc_dim = getfield(structure,'dim');
  elseif isfield(structure,'dimension')
    nc_dim = getfield(structure,'dimension');
  else
    for idim = 1:ndims(data)
      nc_dim(idim) = size(data,idim);    
    end
  end
  % 
  if ~iscell(nc_dim) && isstr(nc_dim) 
    nc_dim = cellstr(nc_dim);
  end 
  % 
  if isnumeric(nc_dim) % new dimensions
    if prod(nc_dim) ~= length(data(:))
      disp('ERROR: assigned dimensions NOT matching data!');
      netcdf.close(nc); 
      status = -1;
      return;
    end
    % def. dimensions
    for idim = 1:length(nc_dim)
      dim_data(idim) = netcdf.defDim(nc,[structure.name,'_dim',num2str(idim)],nc_dim(idim));
    end
  else % existent dimensions
    total_size = 1;
    for idim = 1:length(nc_dim)      
      dim_data(idim) = netcdf.inqDimID(nc,nc_dim{idim});
    end
  end
  % def. variable
  var_data = netcdf.defVar(nc,structure.name,structure.nc_type,[dim_data]);
  % put variable attributes
  if isfield(structure,'attribute') 
    nc_att = getfield(structure,'attribute');
  elseif isfield(structure,'attr') 
    nc_att = getfield(structure,'attr');
  elseif isfield(structure,'att') 
    nc_att = getfield(structure,'att');
  else
    nc_att = [];
  end
  if ~isempty(nc_att)
    att_names = fieldnames(nc_att);
    for iatt = 1:length(att_names)
      netcdf.putAtt(nc,var_data,att_names{iatt},eval(['nc_att.',att_names{iatt}]));
    end
    if isfield(nc_att,'missing_value')
      data(isnan(data)) = nc_att.missing_value;
    end
  end
  % 
  netcdf.endDef(nc);
  % put var
  netcdf.putVar(nc,var_data,data);       
    
elseif strcmpi(structure.type,'dimension') ||  strcmpi(structure.type,'dim')

  % dimension
  if isfield(structure,'if_unlim') && structure.if_unlim == 1
    dim_data = netcdf.defDim(nc,structure.name,netcdf.getConstant('NC_UNLIMITED'));
  else
    dim_data = netcdf.defDim(nc,structure.name,length(data(:)));
  end
  % variable
  var_data = netcdf.defVar(nc,structure.name,structure.nc_type,[dim_data]);
  % put variable attributes
  if isfield(structure,'attribute') 
    nc_att = getfield(structure,'attribute');
  elseif isfield(structure,'attr') 
    nc_att = getfield(structure,'attr');
  elseif isfield(structure,'att') 
    nc_att = getfield(structure,'att');
  else
    nc_att = [];
  end
  if ~isempty(nc_att)
    att_names = fieldnames(nc_att);
    for iatt = 1:length(att_names)
      netcdf.putAtt(nc,var_data,att_names{iatt},eval(['nc_att.',att_names{iatt}]));
    end
    if isfield(nc_att,'missing_value')
      data(isnan(data)) = nc_att.missing_value;
    end
  end
  % 
  netcdf.endDef(nc);
  % put var
  netcdf.putVar(nc,var_data,data);

elseif strcmpi(structure.type,'attribute') ||  strcmpi(structure.type,'att')

  if isfield(structure,'var_name') % variable attribute
    var_data = netcdf.inqVarID(nc,var_name);
    netcdf.putAtt(nc,var_data,structure.name,data);
  else % global attribute
    netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),structure.name,data);
  end  
  netcdf.endDef(nc);

else

  disp('ERROR: Unrecognizable type!');
  status = -1;

end

netcdf.close(nc);
status = 0;

