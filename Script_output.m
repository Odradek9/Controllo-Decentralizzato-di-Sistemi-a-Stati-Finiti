close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                               OUTPUT                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = waitbar(0, "Gathering all sub-stransition systems...");

all_Tc_sub_transition_systems = Split(Tc, Xc0_str, Xcm_str);

all_strategies = cell(length(all_Tc_sub_transition_systems,1), 1);

progress_bar_step = length(all_Tc_sub_transition_systems)/10;
tick = 0;

for atts = 1:length(all_Tc_sub_transition_systems)

    if fix(i/progress_bar_step) > tick
        tick = fix(i/progress_bar_step);
        waitbar((1+(tick*0.1))/2, h ,"Calculating plants stategies for all sub-transition system...");
    end
    
    current_transition_system = all_Tc_sub_transition_systems{atts};

    %----------------------------------------------------------------------
    
    new_Xc = {};
    nc = 1;
    new_Xc_str = strings(0);
    all_cts_states = keys(current_transition_system);
    
    for xc = 1:length(Xc)
    
        xc_state = Xc(xc,:);
        xc_state_str = string(join(cellfun(@(x) num2str(x), xc_state, 'UniformOutput', false), ', '));
    
        if any(strcmp(xc_state_str, all_cts_states))
            new_Xc{nc,1} = xc_state{1};
            new_Xc{nc,2} = xc_state{2};
            new_Xc{nc,3} = xc_state{3};
            new_Xc{nc,4} = xc_state{4};
            nc = nc+1;
    
            new_Xc_str(end+1) = xc_state_str;
        end
    
    end
    
    %----------------------------------------------------------------------
    
    tq1 = specification_transition_tables{1};
    tq2 = specification_transition_tables{2};
    
    size_xc0 = size(Xc0,1);
    
    for xc = 1:size_xc0
    
        xc_fix = size_xc0 - xc + 1;
    
        xc_state = Xc0(xc_fix,:);
        xc_state_str = Xc0_str(xc_fix);
    
        xq1 = xc_state{2};
        values_xq1 = tq1(xq1);
        valid_xq1_succ = values_xq1("successors");
    
        num_of_next = size(valid_xq1_succ,1);
        temp = strings(1, 0);
        tt = 1;
        for nn = 1:num_of_next
            row = valid_xq1_succ{nn,2};
            for mm = 1:length(row)
                temp(tt) = row{mm};
                tt = tt+1;
            end
        end
        valid_xq1_succ = temp;
    
        xq2 = xc_state{4};
        values_xq2 = tq2(xq2);
        valid_xq2_succ = values_xq2("successors");
    
        num_of_next = size(valid_xq2_succ,1);
        temp = strings(1, 0);
        tt = 1;
        for nn = 1:num_of_next
            row = valid_xq2_succ{nn,2};
            for mm = 1:length(row)
                temp(tt) = row{mm};
                tt = tt+1;
            end
        end
        valid_xq2_succ = temp;
    
        values = current_transition_system(xc_state_str);
        succ = values("successors");
    
        size_succ = size(succ,1);
    
        for i = 1:size_succ
            i_fix = size_succ - i + 1;
            next_states = [ succ{i_fix,2} ];
    
            length_next_states = length(next_states);
            
            for j = 1:length_next_states
    
                j_fix = length_next_states - j + 1;
                
                next = next_states{j_fix};
                
                index = find(strcmp(next, new_Xc_str));
                
                if isempty(index)
                    to_remove = true;
                else
    
                    raw_state = new_Xc(index,:);
    
                    to_remove = false;
    
                    if ( ismember(raw_state{2}, valid_xq1_succ) && ismember(raw_state{4}, valid_xq2_succ) )
    
                        values_next_xq1 = tq1(raw_state{2});
                        output_1 = values_next_xq1("output");
    
                        values_next_xq2 = tq2(raw_state{4});
                        output_2 = values_next_xq2("output");
    
                        if ~( isequal(raw_state{1}, output_1) ) ...
                                && ( isequal(raw_state{3}, output_2) ) 
                            to_remove = true;
                        end
    
                    else
                        to_remove = true;
                    end
    
                end
    
                if to_remove
                    
                    next_states(j_fix) = [];
    
                end
    
            end
    
            if isempty(next_states)
                succ(i_fix,:) = [];
            end
    
        end
    
        if isempty(succ)
    
            Xc0(xc_fix,:) = [];
            Xc0_str(xc_fix) = [];
    
        end
    
    end
    
    %----------------------------------------------------------------------
    
    X0_ast = union(intersect(Xc0_str, Xcm_str), Xc0_str);
    
    %----------------------------------------------------------------------
    
    stack = X0_ast;
    checked_states = strings(1,length(new_Xc_str));
    cs = 1;
    finished = 0;
    while ~finished
        if isempty(stack)
            finished = 1;
        elseif ~isKey(current_transition_system, stack(end))
            stack(end) = [];
        else
    
            xc_state_str = stack(end);
    
            checked_states(cs) = stack(end);
            cs = cs+1;
            
            stack(end) = [];
            
            index = find(strcmp(xc_state_str, new_Xc_str));
            xc_state = new_Xc(index,:);
    
            xq1 = xc_state{2};
            values_xq1 = tq1(xq1);
            valid_xq1_succ = values_xq1("successors");
    
            num_of_next = size(valid_xq1_succ,1);
            temp = strings(1, 0);
            tt = 1;
            for nn = 1:num_of_next
                row = valid_xq1_succ{nn,2};
                for mm = 1:length(row)
                    temp(tt) = row{mm};
                    tt = tt+1;
                end
            end
            valid_xq1_succ = temp;
    
            xq2 = xc_state{4};
            values_xq2 = tq2(xq2);
            valid_xq2_succ = values_xq2("successors");
    
            num_of_next = size(valid_xq2_succ,1);
            temp = strings(1, 0);
            tt = 1;
            for nn = 1:num_of_next
                row = valid_xq2_succ{nn,2};
                for mm = 1:length(row)
                    temp(tt) = row{mm};
                    tt = tt+1;
                end
            end
            valid_xq2_succ = temp;
    
            values = current_transition_system(xc_state_str);
            succ = values("successors");
    
            size_succ = size(succ,1);
    
            for i = 1:size_succ
                i_fix = size_succ - i + 1;
                next_states = [ succ{i_fix,2} ];
    
                length_next_states = length(next_states);
    
                for j = 1:length_next_states
    
                    j_fix = length_next_states - j + 1;
    
                    next = next_states{j_fix};
    
                    index = find(strcmp(next, new_Xc_str));
    
                    if isempty(index)
                        to_remove = true;
                    else
    
                        raw_state = new_Xc(index,:);
    
                        to_remove = false;
    
                        if ( ismember(raw_state{2}, valid_xq1_succ) && ismember(raw_state{4}, valid_xq2_succ) ) ...
                                && any(strcmp(next, new_Xc_str))
    
                            values_next_xq1 = tq1(raw_state{2});
                            output_1 = values_next_xq1("output");
    
                            values_next_xq2 = tq2(raw_state{4});
                            output_2 = values_next_xq2("output");
    
                            if ~( isequal(raw_state{1}, output_1) ) ...
                                    && ( isequal(raw_state{3}, output_2) )
                                to_remove = true;
                            else
    
                                if ~any(strcmp(next, checked_states)) && ~any(strcmp(next, stack))
                                    stack(end+1) = next;
                                end
    
                            end
    
                        else
                            to_remove = true;
                        end
    
                    end
    
                    if to_remove
    
                        next_states(j_fix) = [];
    
                    end
    
                end
    
                if isempty(next_states)
                    succ(i_fix,:) = [];
                end
    
            end
    
        end
    end
    
    
    %----------------------------------------------------------------------
    
    
    current_transition_system = Trim(current_transition_system, X0_ast, Xcm_str);
    
    %----------------------------------------------------------------------
    
    Ui_tot = cell(1,2);
    Ui_s_tot = cell(1,2);
    
    for i = 1:2
    
        Ui = containers.Map;
    
        for xc = 1:length(new_Xc)
    
            xc_state = new_Xc(xc,:);
            xc_state_str = new_Xc_str(xc);
    
            all_ui = [];
            ai = 1;
    
            values = current_transition_system(xc_state_str);
            succ = values("successors");
            if ~isempty(succ)
                all_inputs_of_succ = succ(:,1);
                for al = 1:length(all_inputs_of_succ)
        
                    single_input_tuple = all_inputs_of_succ{al};
                    ui = single_input_tuple{i};
        
                    if isempty(all_ui)
                        all_ui(ai, :) = ui;
                        ai = ai+1;
                    elseif ~ismember(ui, all_ui, "rows")
                        all_ui(ai, :) = ui;
                        ai = ai+1;
                    end
        
                end
        
                Ui(xc_state_str) = all_ui;
            end
    
        end
    
        Ui_tot{i} = Ui;
    
    end
    
    %----------------------------------------------------------------------
    
    for i = 1:2
    
        current_Ui = Ui_tot{i};
        
        Ui_s = containers.Map;
    
        for xc = 1:length(new_Xc)
    
            xc_state = new_Xc(xc,:);
            
            xc_state_str = new_Xc_str(xc);
    
            
    
            all_ui_s = [];
            ais = 1;
    
            values = current_transition_system(xc_state_str);
            succ = values("successors");
            for al = 1:size(succ, 1)
    
                all_next_states_for_input = succ{al,2};
    
                for j = 1:length(all_next_states_for_input)
    
                    index = find(strcmp(all_next_states_for_input{j}, new_Xc_str),1);
                    if ~isempty(index)
    
                        str_next = all_next_states_for_input{j};
                        raw_next = new_Xc(index,:);
    
                        if isequal(xc_state{1}, raw_next{1}) && isequal(xc_state{2}, raw_next{2})
    
                            if isKey(current_Ui, str_next)
                                ui_of_plus = current_Ui(str_next);
    
                                for h = 1:length(ui_of_plus)
                                    ui = ui_of_plus(h,:);
    
                                    if isempty(all_ui_s)
                                        all_ui_s(ais, :) = ui;
                                        ais = ais+1;
                                    elseif ~ismember(ui, all_ui_s, "rows")
                                        all_ui_s(ais, :) = ui;
                                        ais = ais+1;
                                    end
    
                                end
    
                            end
    
                        end
    
                    end
    
                end
    
            end
    
            Ui_s(xc_state_str) = all_ui_s;
    
        end
    
        Ui_s_tot{i} = Ui_s;
    
    end
    
    %----------------------------------------------------------------------
    
    U1 = Ui_tot{1};
    U2_s = Ui_s_tot{2};
    
    for i = 1:length(all_cts_states)
    
        state = all_cts_states{i};
        values = current_transition_system(state);
        succ = values("successors");
        pred = values("predecessors");
        
        if isKey(U1, state)
            all_U1_inputs = U1(state);
        else
            all_U1_inputs = [];
        end
    
        if isKey(U2_s, state)
            all_U2_s_inputs = U2_s(state);
        else
            all_U2_s_inputs = [];
        end
        
        if isempty(succ)
            all_inputs_of_succ = [];
        else
            all_inputs_of_succ = succ(:,1);
        end
        size_aios = length(all_inputs_of_succ);
        for al = 1:length(all_inputs_of_succ)
            al_f = size_aios - al + 1;
            single_input_tuple = all_inputs_of_succ{al_f};
    
            if ( isempty(all_U1_inputs) || ~ismember(single_input_tuple{1}, all_U1_inputs, "rows") ) || ...
                    ( isempty(all_U2_s_inputs) || ~ismember(single_input_tuple{2}, all_U2_s_inputs, "rows") )
    
                next_states = succ{al_f,2};
    
                for ns = 1:length(next_states)
                    next = next_states{ns};
    
                    next_values = current_transition_system(next);
                    next_pred = next_values("predecessors");
    
                    input_index = 0;
                    for api = 1:size(next_pred,1)
                        pred_inputs = next_pred{api,1};
                        if isequal(pred_inputs, single_input_tuple)
                            input_index = api;
                            break;
                        end
                    end
    
                    row = next_pred{input_index, 2};
                    str_ar_row = [row{:}];
                    index = find(strcmp(next, str_ar_row));
                    row(index) = [];
                    if isempty(row)
                        next_pred(input_index,:) = [];
                    else
                        next_pred{input_index, 2} = row;
                    end
    
                end
    
                succ(al_f, :) = [];
    
            end
    
        end
    
    end
    
    
    %----------------------------------------------------------------------
    
    current_transition_system = Trim(current_transition_system, Xc0_str, Xcm_str);
    
    
    %----------------------------------------------------------------------
    
    
    all_cts_states = keys(current_transition_system);
    
    Z_tot = cell(1,2);
    
    for p = 1:2
    
        pf = 2*(p-1)+1;
        Zi = cell(0,2);
        z = 1;
    
        plant_i_states = plants_states{p};
        total_specification_states_i = specification_states{p};
    
        for i = 1:size(plant_i_states,1)
            state_i = plant_i_states(i,:);
    
            for qi = 1:length(total_specification_states_i)
                spec_state_i = total_specification_states_i{qi};
    
                tuple = {state_i, spec_state_i};
    
                for ii = 1:length(all_cts_states)
    
                    state = all_cts_states{ii};
                    index = find(strcmp(state, new_Xc_str));
                    raw_state = new_Xc(index,:);
    
                    if ( isequal(raw_state{pf}, tuple{1}) && isequal(raw_state{pf+1}, tuple{2}))
    
                        Zi{z,1} = tuple{1};
                        Zi{z,2} = tuple{2};
                        z = z+1;
    
                        break;
                    end
    
                end
    
            end
    
        end
    
        Z_tot{p} = Zi;
    
    end
    
    
    %----------------------------------------------------------------------
    
    
    Z0_tot = cell(1,2);
    
    for p = 1:2
    
        pf = 2*(p-1)+1;
        Z0i = cell(0,2);
        z = 1;
    
        plant_i_initial_states = plants_initial_states{p};
        total_specification_initial_states_i = specification_initial_states{p};
    
        for i = 1:size(plant_i_initial_states,1)
            state_i = plant_i_initial_states(i,:);
    
            for qi = 1:length(total_specification_initial_states_i)
                spec_state_i = total_specification_initial_states_i{qi};
    
                tuple = {state_i, spec_state_i};
    
                for ii = 1:length(Xc0_str)
    
                    state = Xc0_str{ii};
                    index = find(strcmp(state, Xc0_str));
                    raw_state = Xc0(index,:);
    
                    if ( isequal(raw_state{pf}, tuple{1}) && isequal(raw_state{pf+1}, tuple{2}))
    
                        Z0i{z,1} = tuple{1};
                        Z0i{z,2} = tuple{2};
                        z = z+1;
    
                        break;
                    end
    
                end
    
            end
    
        end
    
        Z0_tot{p} = Z0i;
    
    end
    
    
    %----------------------------------------------------------------------
    
    
    G_tot = cell(0,2);
    
    for p = 1:2
    
        current_Zi = Z_tot{p};
        Gi = {};
        pf = 2*(p-1)+1;
        z = 1;
    
        for i = 1:size(current_Zi, 1)
    
            zi_state = current_Zi(i,:);
    
            all_results = {};
            ar = 1;
    
            for ii = 1:length(all_cts_states)
    
                state_cts = all_cts_states{ii};
                index = find(strcmp(state_cts, new_Xc_str));
                raw_state = new_Xc(index,:);
    
                if ( isequal(raw_state{pf}, zi_state{1}) && isequal(raw_state{pf+1}, zi_state{2}))
    
                    values = current_transition_system(state_cts);
                    succ = values("successors");
    
                    for j = 1:size(succ)
                        all_succ = succ{j,2};
    
                        for k = 1:length(all_succ)
    
                            succ_state = all_succ{k};
                            index_succ = find(strcmp(succ_state, new_Xc_str));
                            raw_succ_state = new_Xc(index_succ,:);
    
                            if isempty(all_results)
                                all_results{ar,1} = raw_succ_state{pf};
                                all_results{ar,2} = raw_succ_state{pf+1};
                                ar = ar+1;
                            else
    
                                new_pair = {raw_succ_state{pf}, raw_succ_state{pf+1}};
    
                                if ~any(cellfun(@(x) isequal(x, new_pair), num2cell(all_results, 2)))
                                    all_results{ar,1} = raw_succ_state{pf};
                                    all_results{ar,2} = raw_succ_state{pf+1};
                                    ar = ar+1;
                                end
                            end
                        end
    
                    end
    
                end
    
            end
    
            condition = {zi_state, zi_state{1}};
    
            Gi{z,1} = condition;
            Gi{z,2} = all_results;
            z = z+1;
    
        end
    
        G_tot{p} = Gi;
    
    end
    
    
    %----------------------------------------------------------------------
    
    
    h_tot = cell(0,2);
    
    for p = 1:2
    
        current_Zi = Z_tot{p};
        hi = {};
        pf = 2*(p-1)+1;
        z = 1;
    
        for i = 1:size(current_Zi, 1)
    
            zi_state = current_Zi(i,:);
    
            all_results = {};
            ar = 1;
    
            for ii = 1:length(all_cts_states)
    
                state_cts = all_cts_states{ii};
                index = find(strcmp(state_cts, new_Xc_str));
                raw_state = new_Xc(index,:);
    
                if ( isequal(raw_state{pf}, zi_state{1}) && isequal(raw_state{pf+1}, zi_state{2}))
    
                    values = current_transition_system(state_cts);
                    succ = values("successors");
    
                    for j = 1:size(succ)
                        inputs = succ{j,1};
    
                        if isempty(all_results)
                            all_results{ar} = inputs{p};
                            ar = ar+1;
                        else
    
                            if ~any(cellfun(@(x) isequal(x, inputs{p}), all_results))
                                all_results{ar} = inputs{p};
                                ar = ar+1;
                            end
                        end
    
                    end
    
                end
    
            end
    
            condition = {zi_state, zi_state{1}};
    
            hi{z,1} = condition;
            hi{z,2} = all_results;
            z = z+1;
    
        end
    
        h_tot{p} = hi;
    
    end
    
    
    %----------------------------------------------------------------------
    
    
    C_tot = cell(1,2);
    
    for p = 1:2
    
        current_Zi = Z_tot{p};
        current_Zi0 = Z0_tot{p};
        current_Gi = G_tot{p};
        current_hi = h_tot{p};
    
       Ci = containers.Map;
    
       for zs = 1:size(current_Zi,1)
    
           zi_state = current_Zi(zs,:);
           zi_state_str =   ;
           values = containers.Map;
           
           if any(cellfun(@(x) isequal(x, zi_state), num2cell(current_Zi0, 2)))
               values("initial") = 1;
           else
               values("initial") = 0;
           end
    
           for hh = 1:size(current_hi,1)
               condition = current_hi{hh,1};
               first_condition = condition{1};
               
               outputs = current_hi{hh,2};
    
               if isequal(first_condition, zi_state)
                   values("output") = outputs;
                    break;
               end
           end
    
           succ_cell_matrix = cell(1,2);
           for kl = 1:size(current_Gi,1)
               condition = current_Gi{kl,1};
               first_condition = condition{1};
               second_condition = condition{2};
               successors = current_Gi{kl,2};
                
               successors_str = strings(1,size(successors,1));
               for l = 1:size(successors,1)
                   row = successors(l,:);
                   row_str = string(join(cellfun(@(x) num2str(x), row, 'UniformOutput', false), ', '));
                   successors_str(l) = row_str;
               end
    
               if isequal(first_condition, zi_state)
                   succ_cell_matrix{1} = second_condition;
                   succ_cell_matrix{2} = successors_str;
                   break;
               end
           end
           values("successors") = succ_cell_matrix;
           Ci(zi_state_str) = values;
       end
    
       %%%
       C_tot{p} = Ci;
    
    end
    
    %----------------------------------------------------------------------
    
    R_tot = cell(1,2);
    
    for p = 1:2
        
        current_Zi0 = Z0_tot{p};
        Ri = cell(size(current_Zi0,1), 1);
    
        for i = 1:size(current_Zi0,1)
    
            zi0_state = current_Zi0(i,:);
            zi0_state_str = string(join(cellfun(@(x) num2str(x), zi0_state, 'UniformOutput', false), ', '));
            Ri{i} = { zi0_state_str, zi0_state{1} };
        
        end
    
        R_tot{p} = Ri;
    
    end
    
    S1 = {C_tot{1}, R_tot{1}};
    S2 = {C_tot{2}, R_tot{2}};
    S_tot = {S1{:}; S2{:}};

    all_strategies{atts} = S_tot;

end

%----------------------------------------------------------------------

waitbar(1, h ,"Job done.");
pause(3);
close(h);

disp("Done.")

clearvars -except all_strategies all_Tc_sub_transition_systems Xc0_str ...
                 Xcm_str Xc
