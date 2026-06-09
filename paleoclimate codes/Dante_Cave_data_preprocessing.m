% 读取Dante Cave石笋数据
filename = 'dantecave2013sletten-noaa.txt'; % 确保文件在当前目录

% % 方法1：使用readtable自动检测格式
% try
%     % 自动检测导入选项
%     opts = detectImportOptions(filename, 'FileType', 'text');
% 
%     % 设置变量名称（根据文件中的列名）
%     opts.VariableNames = {'depth_mm', 'age_calBP', 'calcite_percent', ...
%                          'd13CcarbVPDB', 'd18OcarbVPDB', 'd13CcarbVPDB_ac', ...
%                          'd18OcarbVPDB_ac', 'growth_mmyr', 'notes'};
% 
%     % 设置数据类型
%     opts.VariableTypes = {'double', 'double', 'double', 'double', 'double', ...
%                          'double', 'double', 'double', 'char'};
% 
%     % 设置缺失值处理
%     opts.MissingRule = 'fill';
%     opts.FillValue = NaN;
% 
%     % 跳过注释行（以#开头的行）
%     opts.CommentStyle = {'#'};
% 
%     % 读取数据
%     data = readtable(filename, opts);
% 
% catch ME
%     % 如果自动检测失败，使用手动方法
%     disp('自动检测失败，使用手动读取方法...');
%     data = manualReadData(filename);
% end

data = manualReadData(filename);

% 显示数据基本信息
disp('数据基本信息:');
disp(['数据行数: ', num2str(height(data))]);
disp(['数据列数: ', num2str(width(data))]);
disp(' ');

% 显示列名和数据类型
disp('数据列信息:');
for i = 1:width(data)
    varName = data.Properties.VariableNames{i};
    varType = class(data.(varName));
    
    % 跳过notes列的缺失值检测
    if strcmp(varName, 'notes')
        disp([varName, ' (', varType, ') - 跳过缺失值检测']);
    % 对于其他数值类型列，进行缺失值统计
    elseif isa(data.(varName), 'double')
        nonMissing = sum(~isnan(data.(varName)));
        disp([varName, ' (', varType, ') - 非缺失值: ', num2str(nonMissing), '/', num2str(height(data))]);
    else
        % 对于其他非数值类型，只显示类型信息
        disp([varName, ' (', varType, ') - 跳过缺失值检测']);
    end
end


% 提取前 8 列
data_subset = data(:, 1:8);

% 将表格转换为矩阵
data_sub = table2array(data_subset);

%填补空缺值
filled_data = fill_missing_nearest_mean_matrix(data_sub);




