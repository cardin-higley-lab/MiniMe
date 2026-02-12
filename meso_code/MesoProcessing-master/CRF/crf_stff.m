
roi = 1; 
[crf339_1_loco, se_loco_339_1] = get_crf_responses(crf339_1.data_struct, "loco", roi);
[crf365_1_loco, se_loco_365_1] = get_crf_responses(crf365_1.data_struct, "loco", roi);
[crf365_2_loco, se_loco_365_2] = get_crf_responses(crf365_2.data_struct, "loco", roi);
[crf365_3_loco, se_loco_365_3] = get_crf_responses(crf365_3.data_struct, "loco", roi);
[crf365_4_loco, se_loco_365_4] = get_crf_responses(crf365_4.data_struct, "loco", roi);
control_vals = [crf339_1_loco crf365_1_loco crf365_2_loco  crf365_3_loco crf365_4_loco]; % 
%control_mean = mean(control_vals, 2);
crf_data_v1.control_loco_means = control_vals; 

[crf363_1_loco, se_loco_363_1] = get_crf_responses(crf363_1.data_struct, "loco", roi);
[crf363_2_loco, se_loco_363_2] = get_crf_responses(crf363_2.data_struct, "loco", roi);
[crf366_1_loco, se_loco_366_1] = get_crf_responses(crf366_1.data_struct, "loco", roi);
[crf366_2_loco, se_loco_366_2] = get_crf_responses(crf366_2.data_struct, "loco", roi);
mutant_vals = [crf363_1_loco crf363_2_loco crf366_1_loco crf366_2_loco];
%mutant_mean = mean(mutant_vals, 2);
crf_data_v1.mutant_loco_means = mutant_vals; 


[crf339_1_sit, se_sit_339_1] = get_crf_responses(crf339_1.data_struct, "sit", roi);
[crf365_1_sit, se_sit_365_1] = get_crf_responses(crf365_1.data_struct, "sit", roi);
[crf365_2_sit, se_sit_365_2] = get_crf_responses(crf365_2.data_struct, "sit", roi);
[crf365_3_sit, se_sit_365_3] = get_crf_responses(crf365_3.data_struct, "sit", roi);
[crf365_4_sit, se_sit_365_4] = get_crf_responses(crf365_4.data_struct, "sit", roi);
control_vals = [crf339_1_sit  crf365_1_sit crf365_2_sit crf365_3_sit crf365_4_sit]; %
%control_mean = mean(control_vals, 2);
crf_data_v1.control_sit_means = control_vals; 


[crf363_1_sit, se_sit_363_1] = get_crf_responses(crf363_1.data_struct, "sit", roi);
[crf363_2_sit, se_sit_363_2] = get_crf_responses(crf363_2.data_struct, "sit", roi);
[crf366_1_sit, se_sit_366_1] = get_crf_responses(crf366_1.data_struct, "sit", roi);
[crf366_2_sit, se_sit_366_2] = get_crf_responses(crf366_2.data_struct, "sit", roi);
mutant_vals = [crf363_1_sit crf363_2_sit crf366_1_sit crf366_2_sit];
%mutant_mean = mean(mutant_vals, 2);
crf_data_v1.mutant_sit_means = mutant_vals; 

save("crf_data_v1", "crf_data_v1")




%% modulation index assessment 


for roi = 1:4 
    [crf339_1_loco, se_loco_339_1] = get_crf_responses(crf339_1.data_struct, "loco", roi);
    [crf365_1_loco, se_loco_365_1] = get_crf_responses(crf365_1.data_struct, "loco", roi);
    [crf365_2_loco, se_loco_365_2] = get_crf_responses(crf365_2.data_struct, "loco", roi);
    [crf365_3_loco, se_loco_365_3] = get_crf_responses(crf365_3.data_struct, "loco", roi);
    [crf365_4_loco, se_loco_365_4] = get_crf_responses(crf365_4.data_struct, "loco", roi);
    control_loco_vals = [crf339_1_loco crf365_1_loco crf365_2_loco  crf365_3_loco crf365_4_loco]; % 
    
    [crf339_1_sit, se_sit_339_1] = get_crf_responses(crf339_1.data_struct, "sit", roi);
    [crf365_1_sit, se_sit_365_1] = get_crf_responses(crf365_1.data_struct, "sit", roi);
    [crf365_2_sit, se_sit_365_2] = get_crf_responses(crf365_2.data_struct, "sit", roi);
    [crf365_3_sit, se_sit_365_3] = get_crf_responses(crf365_3.data_struct, "sit", roi);
    [crf365_4_sit, se_sit_365_4] = get_crf_responses(crf365_4.data_struct, "sit", roi);
    control_sit_vals = [crf339_1_sit  crf365_1_sit crf365_2_sit crf365_3_sit crf365_4_sit]; %
    
    [crf363_1_loco, se_loco_363_1] = get_crf_responses(crf363_1.data_struct, "loco", roi);
    [crf363_2_loco, se_loco_363_2] = get_crf_responses(crf363_2.data_struct, "loco", roi);
    [crf366_1_loco, se_loco_366_1] = get_crf_responses(crf366_1.data_struct, "loco", roi);
    [crf366_2_loco, se_loco_366_2] = get_crf_responses(crf366_2.data_struct, "loco", roi);
    mutant_loco_vals = [crf363_1_loco crf363_2_loco crf366_1_loco crf366_2_loco];
    
    [crf363_1_sit, se_sit_363_1] = get_crf_responses(crf363_1.data_struct, "sit", roi);
    [crf363_2_sit, se_sit_363_2] = get_crf_responses(crf363_2.data_struct, "sit", roi);
    [crf366_1_sit, se_sit_366_1] = get_crf_responses(crf366_1.data_struct, "sit", roi);
    [crf366_2_sit, se_sit_366_2] = get_crf_responses(crf366_2.data_struct, "sit", roi);
    mutant_sit_vals = [crf363_1_sit crf363_2_sit crf366_1_sit crf366_2_sit];
    
    control_mod_ind = (control_loco_vals - control_sit_vals)./(control_loco_vals + control_sit_vals);
    mutant_mod_ind = (mutant_loco_vals - mutant_sit_vals)./(mutant_loco_vals + mutant_sit_vals);
    indeces{roi} = [control_mod_ind mutant_mod_ind];
end



figure(); tiledlayout(2,2, "TileSpacing", "tight")
nexttile(); boxchart(indeces{1}); title('V1'); set(gca, "XTickLabel", ["Control", "Mutant"])
nexttile(); boxchart(indeces{2}); title('PM'); set(gca, "XTickLabel", ["Control", "Mutant"])
nexttile(); boxchart(indeces{3}); title('LM'); set(gca, "XTickLabel", ["Control", "Mutant"])
nexttile(); boxchart(indeces{4}); title('M2'); set(gca, "XTickLabel", ["Control", "Mutant"])




