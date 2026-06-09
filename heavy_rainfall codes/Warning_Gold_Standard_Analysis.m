clc;
clear;

%% Determine gold standard from raw data
load('region_data.mat');
load('nanjing.mat');
k=3;
aa='nan';

num_cr = 35;
threshold_percent = 0.05;  
consecutive_threshold = 3; 

data_dims = size(region_data);
num_times1 = data_dims(1);

% CR
data = region_data(45:num_times1, 1, :, :); 

num_times = size(data, 1);
num_points = 16 * 16; 
proportions = zeros(num_times, 1); 

for t = 1:num_times
    current_matrix = squeeze(data(t, 1, :, :)); 

    mask_over = current_matrix > num_cr;
    proportions(t) = sum(mask_over(:)) / num_points;
end

labels = zeros(num_times, 1); 

for t = 1:num_times
    current_matrix = squeeze(data(t, 1, :, :));
    binary_matrix = current_matrix > num_cr;
    [labeled_matrix, num_components] = bwlabel(binary_matrix, 8);
    found_large_component = false;

    for comp = 1:num_components
        component_size = sum(labeled_matrix(:) == comp);

        if component_size >= 10
            found_large_component = true;
            break;
        end
    end

    if found_large_component
        labels(t) = 1;
    else
        labels(t) = 0;
    end
end

% Detect consecutive time periods where proportion >5% for at least 3 points
over_percent = proportions > threshold_percent;
consecutive_over_percent = [];
i = 1;
while i <= num_times
    if over_percent(i)
        start_idx = i;
        count = 0;
        while i <= num_times && over_percent(i)
            count = count + 1;
            i = i + 1;
        end
        if count >= consecutive_threshold
            consecutive_over_percent = [consecutive_over_percent; start_idx];
        end
    else
        i = i + 1;
    end
end

% Detect consecutive time periods where labels == 1 for at least 3 points
consecutive_labels_one = [];
i = 1;
while i <= num_times
    if labels(i) == 1
        start_idx = i;
        count = 0;
        while i <= num_times && labels(i) == 1
            count = count + 1;
            i = i + 1;
        end
        if count >= consecutive_threshold
            consecutive_labels_one = [consecutive_labels_one; start_idx];
        end
    else
        i = i + 1;
    end
end

%% Single-region rainstorm warning analysis
sd_curve = green_y;

if max(sd_curve) <= 60
    fprintf('Note: All raw curve values <=60, no warning signal.\n');

    create_combined_single_axis_plot(sd_curve, [], [], [], proportions, labels, ...
        consecutive_over_percent, consecutive_labels_one, num_cr, threshold_percent, consecutive_threshold,1,k,aa);
    create_combined_single_axis_plot(sd_curve, [], [], [], proportions, labels, ...
        consecutive_over_percent, consecutive_labels_one, num_cr, threshold_percent, consecutive_threshold,0,k,aa);
    
    fprintf('No warning signal');
    
    return;
end

params = struct();
params.smooth_window = 6;         % smoothing window (1 hour)
params.percentile_threshold = 90; 
params.std_multiplier = 2;         
params.plot_results = false;      

[warning_points, threshold, analysis_results] = single_region_warning_combined(sd_curve, params);

fprintf('Combined threshold: %.3f\n', threshold);
fprintf('Number of warning points: %d\n', length(warning_points));

%% 
create_combined_single_axis_plot(sd_curve, threshold, warning_points, analysis_results, ...
    proportions, labels, consecutive_over_percent, consecutive_labels_one, ...
    num_cr, threshold_percent, consecutive_threshold,1,k,aa);
create_combined_single_axis_plot(sd_curve, threshold, warning_points, analysis_results, ...
    proportions, labels, consecutive_over_percent, consecutive_labels_one, ...
    num_cr, threshold_percent, consecutive_threshold,0,k,aa);

%% Create combined single-axis plot
function create_combined_single_axis_plot(sd_curve, threshold, warning_points, analysis_results, ...
    proportions, labels, consecutive_over_percent, consecutive_labels_one, ...
    num_cr, threshold_percent, consecutive_threshold, iii, k, aa)

    figure('Position', [50, 50, 1200, 600], 'Name', 'Combined Rainstorm Warning and Gold-Standard Analysis', ...
           'Color', 'white');

    t1 = 1:length(sd_curve);
    t2 = 1:length(proportions);
    max_time = max(length(sd_curve), length(proportions));

    ax = axes('Position', [0.09, 0.12, 0.86, 0.78]);
    hold on;
    box off;

    upper_min = min(sd_curve);
    upper_max = max(sd_curve);
    upper_range = upper_max - upper_min;

    colors = struct();
    colors.primary_blue = [0.25, 0.58, 0.85];    
    colors.smoothed_curve = [0.65, 0.25, 0.65];  
    colors.threshold_line = [0.20, 0.63, 0.32];  
    colors.warning_point = [0.84, 0.15, 0.16];   
    colors.warning_bg = [0.98, 0.70, 0.68];    

    colors.secondary_blue = [0.58, 0.77, 0.91];  
    colors.primary_orange = [0.95, 0.52, 0.16]; 
    colors.secondary_orange = [0.99, 0.85, 0.60]; 
    colors.primary_green = [0.15, 0.55, 0.34];   
    colors.secondary_green = [0.68, 0.85, 0.71]; 
    colors.purple = [0.20, 0.63, 0.65];         
    colors.gray = [0.50, 0.50, 0.50];           
    colors.light_gray = [0.90, 0.90, 0.90];      
    colors.black = [0.10, 0.10, 0.10];           


    % rainstorm warning curve
    if iii==1
        h1 = plot(t1, sd_curve, '-o', 'Color', colors.primary_blue, 'LineWidth', 4, 'MarkerSize', 8,...
            'DisplayName', 'Raw curve');
    else
        h1 = plot(t1, sd_curve, '-', 'Color', colors.primary_blue, 'LineWidth', 4, ...
            'DisplayName', 'Raw curve');
    end

    % smoothed curve
    if ~isempty(analysis_results)
        smoothed_curve = moving_average_filter_combined(sd_curve, 6);
        h2 = plot(t1, smoothed_curve, '-.', 'Color', colors.smoothed_curve, 'LineWidth', 4.0, ...
            'DisplayName', 'smoothed curve');
    end

    % threshold line
    if ~isempty(threshold)
        h3 = plot([t1(1), t1(end)], [threshold, threshold], '--', ...
            'Color', colors.threshold_line, 'LineWidth', 3.5, ...
            'DisplayName', sprintf('threshold=%.3f', threshold));
    else
        h3 = plot([t1(1), t1(end)], [60, 60], '--', ...
            'Color', colors.secondary_green, 'LineWidth', 3.5, ...
            'DisplayName', 'Reference threshold(60)');
    end

    has_warning_points = ~isempty(warning_points) && any(warning_points);

    if has_warning_points
        h4 = scatter(warning_points, sd_curve(warning_points), ...
            350, colors.warning_point, 'p', 'filled', ...
            'MarkerEdgeColor', [0.3, 0.1, 0.1], 'LineWidth', 1.8, ...
            'DisplayName', 'Warning point');

        for i = 1:min(12, length(warning_points))
            wp = warning_points(i);
            if wp <= length(sd_curve)
                y_value = sd_curve(wp);
                label_y = y_value + upper_range * 0.07;

                text(wp, label_y, sprintf('%d', wp), ...
                    'FontSize', 15, 'FontName', 'Arial', ...
                    'Color', colors.black, 'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom', ...
                    'BackgroundColor', 'none');
            end
        end
    else

        % No rainstorm warning signal
        text_x = max_time * 0.8;
        text_y = upper_max + 0.4 * upper_range; 

        text(text_x, text_y, 'No rainstorm warning signal', ...
            'FontSize', 16, 'Color', colors.warning_point, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'BackgroundColor', [1, 1, 1, 0.7], 'EdgeColor', colors.warning_point, ...
            'LineWidth', 1.5,  ...
            'Margin', 5);
    end

    % Gold-standard marker analysis
    vertical_offset = -1;

    labels_inverted = zeros(size(labels));
    for i = 1:length(labels)
        if labels(i) == 1
            labels_inverted(i) = vertical_offset - 100;
        else
            labels_inverted(i) = vertical_offset;
        end
    end

    h8 = stem(t2, labels_inverted, 'Color', colors.secondary_orange, 'LineWidth', 1.8, ...
        'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'none', ...
        'MarkerEdgeColor', colors.secondary_orange, ...
        'DisplayName', sprintf('>%d adjacent points exceeding threshold', num_cr));

    if ~isempty(consecutive_labels_one)
        marker_y = ones(size(consecutive_labels_one)) * (vertical_offset - 100);
        h9 = scatter(consecutive_labels_one, marker_y, 180, ...
            'MarkerFaceColor', colors.primary_orange, ...
            'MarkerEdgeColor', colors.primary_orange*0.7, ...
            'Marker', 'd', 'LineWidth', 2, ...
            'DisplayName', sprintf('%d consecutive points with label=1 start', consecutive_threshold));

        for i = 1:length(consecutive_labels_one)
            x_pos = consecutive_labels_one(i);
            text(x_pos, marker_y(i) - 8, sprintf('%d', x_pos), ...
                'FontSize', 15, 'FontName', 'Arial', ...
                'Color', colors.primary_orange*0.9, 'FontWeight', 'bold',...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'top', ...
                'BackgroundColor', 'none', 'FontWeight', 'bold');
        end
    end

    proportions_inverted = -proportions * 100 + vertical_offset;

    h5 = plot(t2, proportions_inverted, '-', 'Color', colors.primary_orange, ...
        'LineWidth', 4, 'DisplayName', sprintf('Proportion >%d', num_cr));

    threshold_inverted = -threshold_percent * 100 + vertical_offset;
    h6 = plot([t2(1), t2(end)], [threshold_inverted, threshold_inverted], ':', ...
        'Color', colors.purple, 'LineWidth', 3.5, ...
        'DisplayName', 'Proportion threshold 5%');

    if ~isempty(consecutive_over_percent)
        start_values = proportions_inverted(consecutive_over_percent);
        h7 = scatter(consecutive_over_percent, start_values, 200, ...
            'MarkerFaceColor', colors.purple, ...
            'MarkerEdgeColor', colors.purple*0.7, ...
            'Marker', '^', 'LineWidth', 2, ...
            'DisplayName', sprintf('%d consecutive points >5%% start', consecutive_threshold));

        for i = 1:length(consecutive_over_percent)
            x_pos = consecutive_over_percent(i);
            y_pos = start_values(i);
            text(x_pos, y_pos - 5, sprintf('%d', x_pos), ...
                'FontSize', 15, 'FontName', 'Arial', ...
                'Color', colors.purple*0.9, 'FontWeight', 'bold',...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'top', ...
                'BackgroundColor', 'none', 'FontWeight', 'bold');
        end
    end

    xlabel('time point', 'FontSize', 16, 'FontWeight', 'bold', ...
        'Color', colors.black);
    ylabel('value', 'FontSize', 16, 'FontWeight', 'bold', ...
        'Color', colors.black);

    title('Combined Rainstorm Warning and Gold-Standard Analysis', 'FontSize', 18, 'FontWeight', 'bold', ...
        'Color', colors.black, 'VerticalAlignment', 'bottom');

    y_min_bottom = min([proportions_inverted; labels_inverted]) - 35;
    y_max_top = max(sd_curve) + 80;
    ylim([y_min_bottom, y_max_top]);
    xlim([1, max_time]);

    ax.FontSize = 12;
    ax.LineWidth = 2;
    ax.XColor = colors.black;
    ax.YColor = colors.black;

    legend_items = [];
    legend_labels = {};

    if exist('h1', 'var')
        legend_items = [legend_items, h1];
        legend_labels{end+1} = get(h1, 'DisplayName');
    end

    if exist('h2', 'var')
        legend_items = [legend_items, h2];
        legend_labels{end+1} = get(h2, 'DisplayName');
    end

    if exist('h3', 'var')
        legend_items = [legend_items, h3];
        legend_labels{end+1} = get(h3, 'DisplayName');
    end

    if has_warning_points && exist('h4', 'var')
        legend_items = [legend_items, h4];
        legend_labels{end+1} = 'warning point';
    end

    if exist('h5', 'var')
        legend_items = [legend_items, h5];
        legend_labels{end+1} = get(h5, 'DisplayName');
    end

    if exist('h6', 'var')
        legend_items = [legend_items, h6];
        legend_labels{end+1} = get(h6, 'DisplayName');
    end

    if exist('h7', 'var')
        legend_items = [legend_items, h7];
        legend_labels{end+1} = get(h7, 'DisplayName');
    end

    if exist('h8', 'var')
        legend_items = [legend_items, h8];
        legend_labels{end+1} = get(h8, 'DisplayName');
    end

    if exist('h9', 'var')
        legend_items = [legend_items, h9];
        legend_labels{end+1} = get(h9, 'DisplayName');
    end

    plot([1, max_time], [vertical_offset - 80, vertical_offset - 80], '--', ...
        'Color', colors.gray, 'LineWidth', 1, 'HandleVisibility', 'off');


    if has_warning_points
        for i = 1:length(warning_points)
            wp = warning_points(i);
            if wp <= max_time

                h_ref = plot([wp, wp], [y_min_bottom, y_max_top], ...
                    ':', 'Color', [0.8, 0.5, 0.5], 'LineWidth', 3.5);

                uistack(h_ref, 'bottom');

                if i == 1
                    set(h_ref, 'DisplayName', 'Warning point vertical reference');
                else
                    set(h_ref, 'HandleVisibility', 'off');
                end
            end
        end
    end

    ax.Color = 'white';

    ax.Box = 'off';
    ax.Layer = 'top';

    set(gcf, 'PaperPositionMode', 'auto');
    set(gcf, 'InvertHardcopy', 'off');

    hold off;
end

function [warning_points, threshold, analysis_results] = single_region_warning_combined(sd_curve, params)


    n_points = length(sd_curve);


    if any(isnan(sd_curve))
        error('Data contains NaN values; please handle missing data first');
    end

    raw_max = max(sd_curve);

    smoothed_curve = moving_average_filter_combined(sd_curve, params.smooth_window);

    basic_stats.mean = mean(sd_curve);
    basic_stats.median = median(sd_curve);
    basic_stats.std = std(sd_curve);
    basic_stats.min = min(sd_curve);
    basic_stats.max = max(sd_curve);

    %% multiple thresholds
    smoothed_stats.mean = mean(smoothed_curve);
    smoothed_stats.median = median(smoothed_curve);
    smoothed_stats.std = std(smoothed_curve);
    n = length(smoothed_curve);

    % Percentile
    threshold1 = prctile(smoothed_curve, params.percentile_threshold);

    % Mean + standard deviation multiplier
    threshold2 = smoothed_stats.mean + params.std_multiplier * smoothed_stats.std;

    % Median + MAD
    mad_value = median(abs(smoothed_curve - smoothed_stats.median));
    if mad_value < eps
        threshold3 = smoothed_stats.median;
    else
        threshold3 = smoothed_stats.median + 1.5 * 1.4826 * mad_value;
    end

    % Sliding window maximum
    window_size = min(12, floor(n/10));
    if window_size < 3
        window_size = 3;
    end
    window_max_vals = zeros(1, n - window_size + 1);
    for i = 1:(n - window_size + 1)
        window_max_vals(i) = max(smoothed_curve(i:i+window_size-1));
    end
    threshold4 = mean(window_max_vals);

    % Based on diurnal variation
    points_per_day = 144;
    if n >= points_per_day
        num_full_days = floor(n / points_per_day);
        daily_means = zeros(1, num_full_days);
        for d = 1:num_full_days
            start_idx = (d-1)*points_per_day + 1;
            end_idx = d*points_per_day;
            daily_means(d) = mean(smoothed_curve(start_idx:end_idx));
        end
        baseline = mean(daily_means);
        amplitude = (max(daily_means) - min(daily_means)) / 2;
        threshold5 = baseline + 2 * amplitude;
    else
        threshold5 = mean([threshold1, threshold2, threshold3, threshold4]);
    end

    % EWMA control chart
    lambda = 0.2;
    ewma_series = zeros(1, n);
    ewma_series(1) = smoothed_curve(1);
    for i = 2:n
        ewma_series(i) = lambda * smoothed_curve(i) + (1-lambda) * ewma_series(i-1);
    end
    residual_std = std(smoothed_curve - ewma_series);
    threshold6 = mean(ewma_series) + 2 * residual_std;

    %% Combined threshold
    threshold_values = [threshold1, threshold2, threshold3, threshold4, threshold5, threshold6];
    valid_idx = ~isnan(threshold_values);
    valid_thresholds = threshold_values(valid_idx);

    weights = [0.30, 0.15, 0.05, 0.25, 0.05, 0.20];
    weights = weights(valid_idx);
    weights = weights / sum(weights);

    weighted_threshold = sum(valid_thresholds .* weights);
    final_threshold = weighted_threshold;

    %% Threshold rationality correction
    golden_standard = 60;

    if raw_max > golden_standard
        prct95 = prctile(sd_curve, 98);

        if final_threshold < golden_standard
            threshold = max(golden_standard, prct95);
        else
            threshold = final_threshold;
        end
    else
        threshold = final_threshold;
    end

    %% Detect warning points
    above_threshold = smoothed_curve > threshold;

    logical_vec = above_threshold(:);
    padded = [0; logical_vec; 0];
    diff_vec = diff(padded);
    start_idx = find(diff_vec == 1);

    warning_points = [];
    if ~isempty(start_idx)
        warning_points = start_idx';
    end

    analysis_results = struct();
    analysis_results.threshold = threshold;
    analysis_results.raw_max = raw_max;
    analysis_results.golden_standard = golden_standard;
    analysis_results.basic_stats = basic_stats;
    analysis_results.warning_points = warning_points;
    analysis_results.has_warning_points = ~isempty(warning_points); 
end

%% 
function smoothed = moving_average_filter_combined(data, window_size)
    kernel = ones(1, window_size) / window_size;
    smoothed = conv(data, kernel, 'same');
    half_win = floor(window_size/2);
    smoothed(1:half_win) = data(1:half_win);
    smoothed(end-half_win+1:end) = data(end-half_win+1:end);
end

%% 
function [storm_events, event_durations, event_intensities] = identify_storm_events(data, data_type, merge_interval, threshold_percent)
   
    if strcmp(data_type, 'A')
        binary_data = data(:);
    else
        binary_data = data > threshold_percent;
    end
  
    diff_data = diff([0; binary_data; 0]);
    start_idx = find(diff_data == 1);
    end_idx = find(diff_data == -1) - 1;

    events = [start_idx, end_idx];

    if size(events, 1) > 1
        merged_events = [];
        current_event = events(1, :);
        
        for i = 2:size(events, 1)
            if (events(i, 1) - current_event(2)) <= merge_interval
                current_event(2) = events(i, 2);
            else
                merged_events = [merged_events; current_event];
                current_event = events(i, :);
            end
        end
        merged_events = [merged_events; current_event];
        events = merged_events;
    end

    num_events = size(events, 1);
    event_durations = zeros(num_events, 1);
    event_intensities = zeros(num_events, 1);
    
    for i = 1:num_events
        event_durations(i) = events(i, 2) - events(i, 1) + 1;
        
        if strcmp(data_type, 'A')

            event_intensities(i) = 1;
        else

            event_intensities(i) = mean(data(events(i, 1):events(i, 2)));
        end
    end
    
    storm_events = events;
end

