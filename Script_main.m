clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                 MAIN                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------

%%% CREAZIONE STATI Xc

h = waitbar(0, "Calculating all valid Xc states...");

Xc = {};
Xc_str = strings(0);
c = 1;

plant_1_states = plants_states{1};
plant_2_states = plants_states{2};

total_specification_states_1 = specification_states{1};
total_specification_states_2 = specification_states{2};

specification_1_output_matrix = specification_state_output_cell_matrices{1};
specification_2_output_matrix = specification_state_output_cell_matrices{2};

for i = 1:size(plant_1_states,1)
    state_1 = plant_1_states(i,:);

    for j = 1:length(plant_2_states)
        state_2 = plant_2_states(j,:);

        for q1 = 1:length(total_specification_states_1)
            spec_state_1 = total_specification_states_1{q1};
            index_1 = find([specification_1_output_matrix{:,1}] == spec_state_1);
            spec_state_1_output = specification_1_output_matrix{index_1, 2};

            for q2 = 1:length(total_specification_states_2)
                spec_state_2 = total_specification_states_2{q2};
                index_2 = find([specification_2_output_matrix{:,1}] == spec_state_2);
                spec_state_2_output = specification_2_output_matrix{index_2, 2};
                
                if ( d(state_1, spec_state_1_output) <= tolerances(1) ) && ( d(state_2, spec_state_2_output) <= tolerances(2) )
                    Xc{c,1} = state_1;
                    Xc{c,2} = spec_state_1;
                    Xc{c,3} = state_2;
                    Xc{c,4} = spec_state_2;
                    temp = {state_1, spec_state_1, state_2, spec_state_2};
                    Xc_str(c) = string(join(cellfun(@(x) num2str(x), temp, 'UniformOutput', false), ', '));
                    c = c+1;
                end
            end

        end

    end

end

%--------------------------------------------------------------------------

%%% CREAZIONE STATI INIZIALI Xc0

waitbar(1/8, h ,"Calculating all Xc initial state...");

plant_1_initial_states = plants_initial_states{1};
plant_2_initial_states = plants_initial_states{2};
spec_initial_states_1 = specification_initial_states{1};
spec_initial_states_2 = specification_initial_states{2};

Xc0 = {};
Xc0_str = strings(0);
c = 1;

for i = 1:size(plant_1_initial_states,1)

    for j = 1:size(plant_2_initial_states,1)

        for q1_0 = 1:length(spec_initial_states_1)

            for q2_0 = 1:length(spec_initial_states_2)
                
                temp = {plant_1_initial_states(i,:), spec_initial_states_1(q1_0), plant_2_initial_states(j,:), spec_initial_states_2(q2_0)};
                Xc0_index = 0;
                for ii = 1:size(Xc, 1)
                    if isequal(Xc(ii, :), temp)
                        Xc0_index = 1;
                        break;
                    end
                end
                if any(Xc0_index)
                    Xc0{c,1} = temp{1};
                    Xc0{c,2} = temp{2};
                    Xc0{c,3} = temp{3};
                    Xc0{c,4} = temp{4};
                    Xc0_str(c) = string(join(cellfun(@(x) num2str(x), temp, 'UniformOutput', false), ', '));
                    c = c+1;
                end
            end

        end

    end

end

%--------------------------------------------------------------------------

%%% CREAZIONE INPUT Uc

waitbar(2/8, h ,"Calculating all inputs Uc...");

inputs_1 = plants_inputs{1};
inputs_2 = plants_inputs{2};
Uc = {};
c = 1;
for i = 1:length(inputs_1)
    for j = 1:length(inputs_2)
        Uc{c,1} = inputs_1(i,:);
        Uc{c,2} = inputs_2(j,:);
        c = c+1;
    end
end

%--------------------------------------------------------------------------

%%% CREAZIONE STATI MARCATI Xm0

waitbar(3/8, h ,"Calculating all Xc marked states...");

Xcm = {};
Xcm_str = strings(0);
c = 1;
specification_1_marked_states = specification_marked_states{1};
specification_2_marked_states = specification_marked_states{2};

for i = 1:size(Xc,1)
    tuple = { Xc{i,:} };    
    if (any(strcmp(specification_1_marked_states, tuple{2}))) && (any(strcmp(specification_2_marked_states, tuple{4})))
        Xcm{c,1} = tuple{1};
        Xcm{c,2} = tuple{2};
        Xcm{c,3} = tuple{3};
        Xcm{c,4} = tuple{4};
        Xcm_str(c) = string(join(cellfun(@(x) num2str(x), tuple, 'UniformOutput', false), ', '));
        c = c+1;
    end
end

%--------------------------------------------------------------------------

%%% CREAZIONE USCITE Yc

waitbar(4/8, h ,"Calculating all outputs Yc...");

Yc = {};
c = 1;
p1_states = plants_states{1};
p2_states = plants_states{2};
for i = 1:size(p1_states, 1)
    for j = 1:size(p2_states, 1)
        Yc{c,1} = p1_states(i,:);
        Yc{c,2} = p2_states(j,:);
        c = c+1;
    end
end

%--------------------------------------------------------------------------

%%% CREAZIONE CELL MATRICE FUNZIONE DI OUTPUT Hc

waitbar(5/8, h ,"Calculating state-output function Hc...");

Hc = cell(size(Xc,1), 2);
c = 1;
Hq1 = specification_state_output_cell_matrices{1};
Hq2 = specification_state_output_cell_matrices{2};
for i = 1:size(Xc,1)
    tuple = { Xc{i,:} };
    Hc{i,1} = tuple;

    % Parentesi quadre per il find
    Hq1_spec_states = [Hq1{:,1}];

    index_1 = find(Hq1_spec_states == tuple{2});

    % Parentesi quadre per il find
    Hq2_spec_states = [Hq2{:,1}];

    index_2 = find(Hq1_spec_states == tuple{4});
    Hc{i,2} = {Hq1{index_1, 2}, Hq2{index_2, 2}};
end

%--------------------------------------------------------------------------

%%% CREAZIONE SISTEMA DI TRANSIZIONE Tc

waitbar(6/8, h ,"Processing transition system Tc...");

P1_t_s = plants_transition_tables{1};
P2_t_s = plants_transition_tables{2};
Q1_t_s = specification_transition_tables{1};
Q2_t_s = specification_transition_tables{2};

Tc = containers.Map;

progress_bar_step = size(Xc,1)/10;
tick = 0;

for i = 1:size(Xc,1)
    
    if fix(i/progress_bar_step) > tick
        tick = fix(i/progress_bar_step);
        waitbar((6+(tick*0.1))/8, h ,"Processing transition system Tc...");
    end
    
    Xc_state = { Xc{i,:} };
    Xc_state_str = string(join(cellfun(@(x) num2str(x), Xc_state, 'UniformOutput', false), ', '));

    if isKey(Tc, Xc_state_str)
        c_values = Tc(Xc_state_str);
        Xc_state_pred = c_values("predecessors");
    else
        c_values = containers.Map;
        Xc_state_pred = cell(0, 2);
    end
    
    output = Hc(i, 2);
    c_values("output") = output;
    
    c_values("marked") = 0;
    for j = 1:length(Xcm)
        if isequal(Xcm{j}, Xc_state)
            c_values("marked") = 1;
            break
        end
    end
    
    c_values("initial") = 0;
    for j = 1:length(Xc0)
        if isequal(Xc0{j}, Xc_state)
            c_values("initial") = 1;
            break
        end
    end

    x1_state = Xc_state{1};
    x2_state = Xc_state{3};
    x1_spec_state = Xc_state{2};
    x2_spec_state = Xc_state{4};

    x1_values = P1_t_s(num2str(x1_state));
    x1_succ = x1_values("successors");

    x2_values = P2_t_s(num2str(x2_state));
    x2_succ = x2_values("successors");

    x1_spec_values = Q1_t_s(x1_spec_state);
    x1_spec_succ = x1_spec_values("successors");
    if ~isempty(x1_spec_succ)
        x1_spec_all_succ = x1_spec_succ{:,2};
        x1_spec_all_succ = [ x1_spec_all_succ{:,:} ];
        x1_spec_succ = unique(x1_spec_all_succ);
    end

    x2_spec_values = Q2_t_s(x2_spec_state);
    x2_spec_succ = x2_spec_values("successors");
    if ~isempty(x2_spec_succ)
        x2_spec_all_succ = x2_spec_succ{:,2};
        x2_spec_all_succ = [ x2_spec_all_succ{:,:} ];
        x2_spec_succ = unique(x2_spec_all_succ);
    end
    
    Xc_state_succ = {};
    c = 1;

    for u = 1:size(Uc,1)

        u1 = Uc{u,1};
        u2 = Uc{u,2};

        x1_cut_succ = [];
    
        x1_check_condition = {x1_state, x2_state, u1};
        if isempty(x1_succ)
            x1_index = 0;
        else
            all_x1_succ_conditions = [ x1_succ(:,1) ];
            x1_index = 0;
            for g = 1:size(all_x1_succ_conditions, 1)
                if isequal(all_x1_succ_conditions{g, :}, x1_check_condition)
                    x1_index = g;
                    break;
                end
            end
        end

        if x1_index ~= 0

            temp = x1_succ{x1_index,2};
            x1_cut_succ = cell(size(temp, 1),1);
            for xs = 1:size(temp, 1)
                x1_cut_succ{xs} = temp(xs,:);
            end
        end

        x2_cut_succ = [];
    
        x2_check_condition = {x2_state, x1_state, u2};
        if isempty(x2_succ)
            x2_index = 0;
        else
            all_x2_succ_conditions = [ x2_succ(:,1) ];
            x2_index = 0;
            for g = 1:size(all_x2_succ_conditions, 1)
                if isequal(all_x2_succ_conditions{g, :}, x2_check_condition)
                    x2_index = g;
                    break;
                end
            end
        end
        
        if x2_index ~= 0
            
            temp = x2_succ{x2_index,2};
            x2_cut_succ = cell(size(temp, 1),1);
            for xs = 1:size(temp, 1)
                x2_cut_succ{xs} = temp(xs,:);
            end
        end

        Xc_plus_total = {};

        for x1 = 1:size(x1_cut_succ,1)
            x1_plus = x1_cut_succ(x1,:);

            for x2 = 1:size(x2_cut_succ,1)
                x2_plus = x2_cut_succ(x2,:);
                
                for x1q = 1:length(x1_spec_succ)
                    x1q_plus = x1_spec_succ(x1q);
                    
                    for x2q = 1:length(x2_spec_succ)
                        x2q_plus = x2_spec_succ(x2q);

                        Xc_plus = {x1_plus{1}, x1q_plus, x2_plus{1}, x2q_plus};
                        
                        Xc_plus_str = string(join(cellfun(@(x) num2str(x), Xc_plus, 'UniformOutput', false), ', '));
                        Xc_plus_index = strcmp(Xc_plus_str, Xc_str);

                        if any(Xc_plus_index)

                            Xc_plus_total{end+1} = Xc_plus_str;

                            if isequal(Xc_plus_str, Xc_state_str)

                                if isempty(Xc_state_pred)

                                    Xc_state_pred{1,1} = {u1, u2};
                                    Xc_state_pred{1,2} = {Xc_state_str};

                                else

                                    input_index = 0;
                                    for api = 1:size(Xc_state_pred,1)
                                        pred_inputs = Xc_state_pred{api,1};
                                        if isequal(pred_inputs, {u1, u2})
                                            input_index = api;
                                            break;
                                        end
                                    end
                                    if input_index == 0
                                        next_index = size(Xc_state_pred,1) + 1;
                                        Xc_state_pred{next_index,1} = {u1, u2};
                                        Xc_state_pred{next_index,2} = {Xc_state_str};
                                    else
                                        row = Xc_state_pred{input_index, 2};
                                        row{end+1} = Xc_state_str;
                                        Xc_state_pred{input_index, 2} = row;
                                    end

                                end

                            else

                                if isKey(Tc, Xc_plus_str)

                                    plus_values = Tc(Xc_plus_str);
                                    Xc_plus_pred = plus_values("predecessors");
                                    if isempty(Xc_plus_pred)
                                        Xc_plus_pred{1,1} = {u1, u2};
                                        Xc_plus_pred{1,2} = {Xc_state_str};
                                    else

                                        input_index = 0;
                                        for api = 1:size(Xc_plus_pred,1)
                                            pred_inputs = Xc_plus_pred{api,1};
                                            if isequal(pred_inputs, {u1, u2})
                                                input_index = api;
                                                break;
                                            end
                                        end
                                        if input_index == 0
                                            next_index = size(Xc_plus_pred,1) + 1;
                                            Xc_plus_pred{next_index,1} = {u1, u2};
                                            Xc_plus_pred{next_index,2} = {Xc_state_str};
                                        else
                                            row = Xc_plus_pred{input_index, 2};
                                            row{end+1} = Xc_state_str;
                                            Xc_plus_pred{input_index, 2} = row;
                                        end

                                    end

                                    plus_values("predecessors") = Xc_plus_pred;

                                else
                                    plus_values = containers.Map;
                                    Xc_plus_pred = cell(1,2);
                                    Xc_plus_pred{1,1} = {u1, u2};
                                    Xc_plus_pred{1,2} = {Xc_state_str};
                                    plus_values("predecessors") = Xc_plus_pred;

                                    Tc(Xc_plus_str) = plus_values;
                                end

                            end

                            
                        end

                    end

                end

            end

        end
        if ~isempty(Xc_plus_total)
            Xc_state_succ{c,1} = {u1, u2};
            Xc_state_succ{c,2} = Xc_plus_total;
            c = c+1;
        end

    end
    
    c_values("successors") = Xc_state_succ;
    c_values("predecessors") = Xc_state_pred;

    Tc(Xc_state_str) = c_values;

end

waitbar(7/8, h ,"Trimming transition system Tc...");

Tc = Trim(Tc, Xc0_str, Xcm_str);

waitbar(1, h ,"Job done.");
pause(3);
close(h);

disp("Done.")

clearvars -except Tc Xc Xc0 Hc Yc Xcm Uc Xc0_str Xcm_str Xc_str ...
    plants_states plants_initial_states plants_inputs ...
    plants_state_transitions plants_transition_tables ...
    specification_initial_states specification_marked_states ...
    specification_transition_tables specification_transition_matrices ...
    specification_state_output_cell_matrices tolerances ...
    specification_inputs specification_states

%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              FUNZIONI                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Funzione metrica d
function distance = d(m, n)

    distance = norm(m-n);

end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~