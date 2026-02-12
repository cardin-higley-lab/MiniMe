function [mean_responses, se_responses] = get_loco_v1responses(struct, state)

v1_responses = struct.v1_responses;
loco = struct.loco;

switch state
    case "loco"
        v1_responses_loco = nan(size(v1_responses));
        v1_responses_loco(:, loco==1) = v1_responses(:, loco==1);
        mean_responses = mean(v1_responses_loco, 2, 'omitnan');
        se_responses = std(v1_responses_loco, 0, 2, 'omitnan')/sqrt(size(v1_responses_loco, 2));
    case "sit"
        v1_responses_loco = nan(size(v1_responses));
        v1_responses_loco(:, loco==0) = v1_responses(:, loco==0);
        mean_responses = mean(v1_responses_loco, 2, 'omitnan');
        se_responses = std(v1_responses_loco, 0, 2, 'omitnan')/sqrt(size(v1_responses_loco, 2));
end