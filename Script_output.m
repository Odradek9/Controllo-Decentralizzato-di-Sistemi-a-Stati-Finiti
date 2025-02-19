close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                               OUTPUT                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;

%all_Tc_sub_transition_systems = Split(Tc, Xc0_str, Xcm_str)
%Split(specification_transition_tables{1}, specification_initial_states{1}, specification_marked_states{1})

current_transition_system = Tc;

%X0_ast = union(intersect(Xc0_str, Xcm_str), Xc0_str);

new_Xc = {};
nc = 1;
new_Xc_str = strings(0);
all_cts_states = keys(current_transition_system);
in_index = find(strcmp("inputs", all_cts_states));
all_cts_states(in_index) = [];

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

Ui_tot = cell(1,2);
Ui_s_tot = cell(1,2);

for i = 1:2

    %Ui = cell(0,2);
    Ui = containers.Map;

    for xc = 1:length(new_Xc)

        xc_state = new_Xc(xc,:);
        xc_state_str = string(join(cellfun(@(x) num2str(x), xc_state, 'UniformOutput', false), ', '));

        %Ui{xc,1} = xc_state_str;

        all_ui = [];
        ai = 1;

        values = current_transition_system(xc_state_str);
        succ = values("successors");
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

        %Ui{xc,2} = all_ui;
        Ui(xc_state_str) = all_ui;

    end

    %disp(Ui)
    Ui_tot{i} = Ui;

end

for i = 1:2

    current_Ui = Ui_tot{i};
    %Ui_s = cell(0,2);
    Ui_s = containers.Map;

    for xc = 1:length(new_Xc)

        xc_state = new_Xc(xc,:);
        xc_state_str = string(join(cellfun(@(x) num2str(x), xc_state, 'UniformOutput', false), ', '));

        %Ui_s{xc,1} = xc_state_str;

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

        %Ui_s{xc,2}  = all_ui_s;
        Ui_s(xc_state_str) = all_ui_s;

    end

    Ui_s_tot{i} = Ui_s;

end

toc;
tic;

U1 = Ui_tot{1};
U2_s = Ui_s_tot{2};

for i = 1:length(all_cts_states)

    state = all_cts_states{i};
    values = current_transition_system(state);
    succ = values("successors");
    pred = values("predecessors");

    all_U1_inputs = U1(state);
    all_U2_s_inputs = U2_s(state);

    all_inputs_of_succ = succ(:,1);
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

                %all_pred_inputs = next_pred{:,1};
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

toc;
tic;

current_transition_system = Trim(current_transition_system, Xc0_str, Xcm_str);

toc;
tic;

all_cts_states = keys(current_transition_system);
in_index = find(strcmp("inputs", all_cts_states));
all_cts_states(in_index) = [];

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

toc;
tic;

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

toc;
tic;

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
                            %%% DA RIFARE

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

        Gi{z,1} = zi_state;
        Gi{z,2} = zi_state{1};
        Gi{z,3} = all_results;
        z = z+1;

    end

    G_tot{p} = Gi;

end

toc;
tic;

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

        hi{z,1} = zi_state;
        hi{z,2} = zi_state{1};
        hi{z,3} = all_results;
        z = z+1;

    end

    h_tot{p} = hi;

end

toc;
tic;

C_tot = cell(1,2);

for p = 1:2

    current_Zi = Z_tot{p};
    current_Zi0 = Z0_tot{p};
    current_Gi = G_tot{p};
    current_hi = h_tot{p};

   Ci = containers.Map;
   %pi_transition_table("inputs") = pi_inputs;

   for zs = 1:size(current_Zi,1)

       zi_state = current_Zi(zs,:);
       zi_state_str = string(join(cellfun(@(x) num2str(x), zi_state, 'UniformOutput', false), ', '));
       values = containers.Map;
       
       if any(cellfun(@(x) isequal(x, zi_state), num2cell(current_Zi0, 2)))
           values("initial") = 1;
       else
           values("initial") = 0;
       end

       for hh = 1:size(current_hi,1)
           first_condition = current_hi{hh,1};
           second_condition = current_hi{hh,2};
           outputs = current_hi{hh,3};

           if isequal(first_condition, zi_state)
               values("output") = outputs;
                break;
           end
       end

       succ_cell_matrix = cell(1,2);
       for kl = 1:size(current_Gi,1)
           first_condition = current_Gi{kl,1};
           second_condition = current_Gi{kl,2};
           successors = current_Gi{kl,3};

           if isequal(first_condition, zi_state)
               succ_cell_matrix{1} = second_condition;
               succ_cell_matrix{2} = successors;
               break;
           end
       end
       values("successors") = succ_cell_matrix;
       Ci(zi_state_str) = values;
   end

   %%%
   C_tot{p} = Ci;

end

toc;
tic;

R_tot = cell(1,2);

for p = 1:2
    
    current_Zi0 = Z0_tot{p};
    Ri = cell(size(current_Zi0,1), 2);

    for i = 1:size(current_Zi0,1)

        zi0_state = current_Zi0(i,:);
        zi0_state_str = string(join(cellfun(@(x) num2str(x), zi0_state, 'UniformOutput', false), ', '));
        Ri{i,1} = zi0_state_str;
        Ri{i,2} = zi0_state{1};
    
    end

    R_tot{p} = Ri;

end

toc;
tic;

S1 = {C_tot{1}, R_tot{1}};
S2 = {C_tot{2}, R_tot{2}};
S_tot = {S1, S2};

toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              FUNZIONI                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Funzione per ottenere tutti i sotto-sistemi di transizione di un dato
% transition system
function result = Split(transition_table, initial_states, marked_states)

%k = keys(transition_table);

%%% RICERCA DEI PATH

% Andiamo come prima cosa ad implementare una ricerca DFS modificata,
% dove ci segniamo i tragitti che finiscono su uno stato marcato, o che
% contengono uno stato marcato e finiscono su se stessi (in un punto
% prima dello stato marcato)
all_paths = {};
for in = 1:length(initial_states)
    for mar = 1:length(marked_states)

        % Inizializzazione ricerca
        stack = initial_states(in);
        path = strings(0);
        all_visited_states = strings(0);
        coac_states = strings(0);
        index_of_marked = 0;
        branches = 0;
        last_branch_index = 0;
        visited_marked = 0;

        while true

            if isempty(stack)
                break;
            end

            path(end+1) = stack(end);

            % Se finiamo sullo stato marcato ci segnamo il path ad il
            % suo indice su questo path
            if isequal(path(end), marked_states(mar))
                index_of_marked = length(path);
                all_paths{end+1} = path;
            end

            value = transition_table( stack(end) );
            successors = value("successors");

            num_of_next = size(successors,1);
            temp = strings(1, 0);
            tt = 1;
            for nn = 1:num_of_next
                row = successors{nn,2};
                for mm = 1:length(row)
                    temp(tt) = row{mm};
                    tt = tt+1;
                end
            end
            possible_next_states = temp;

            % Non aggiungiamo lo stato marcato alla lista degli stati
            % visitati per permetterci di continuare i percorsi oltre
            % la visita dello stato marcato
            if ~isequal(stack(end), marked_states(mar))
                all_visited_states(end+1) = stack(end);
            end

            temp_coac = coac_states;
            temp_visited = all_visited_states;
            % Se abbiamo già visitato lo stato marcato una volta, lo
            % includiamo nel controllo delle ramificazioni.
            % Non vogliamo contare transizioni su di esso come una
            % transizione valida se ci siamo già stati, altrimenti
            % quando dobbiamo retrocedere su una ramificazione
            % precedente avremmo problemi
            if visited_marked
                temp_coac(end+1) = marked_states(mar);
                temp_visited(end+1) = marked_states(mar);
            end
            A = intersect(possible_next_states, temp_coac);

            % Controllo e preservazione dell'indice, all'interno del
            % path attuale, delle ramificazioni presenti (non contando
            % transizioni su se stesso e su stati antecedenti allo
            % stato marcato)
            self_loop_check = any(strcmp(stack(end), possible_next_states));
            if length(possible_next_states) - self_loop_check - length(A) > 1
                new_index = find(strcmp( stack(end), path), 1);
                if ~ismember(new_index, branches)
                    branches(end+1) = find(strcmp( stack(end), path), 1);
                end
            end

            if isempty(branches)
                last_branch_index = 0;
            else
                last_branch_index = branches(end);
            end

            stack(end) = [];
            if ~isempty(possible_next_states) && ~(all(ismember(possible_next_states, temp_visited)))

                for i = 1:length(possible_next_states)

                    if ~any(strcmp(possible_next_states(i), all_visited_states))

                        % Se ancora non abbiamo incontrato lo stato
                        % marcato in questo path lo possiamo aggiungere
                        % allo stack così da andare a controllare gli
                        % stati successivi oltre di esso
                        if ~(isequal(possible_next_states(i), marked_states(mar)) && visited_marked)
                            stack(end+1) = possible_next_states(i);
                        end

                        % Se il prossimo stato di riferimento è uno
                        % stato antecedente allo stato marcato e
                        % abbiamo già visitato lo stato marcato una
                        % volta allora ci segnamo il percorso attuale e
                        % aggiungiamo gli stati nuovi del path agli
                        % stati antecedenti al marcato (marcato non
                        % incluso)
                        if any(strcmp(possible_next_states(i), coac_states)) && visited_marked
                            all_paths{end+1} = path;
                            coac_states = union( coac_states, path( (index_of_marked+1):end ) );
                        end

                        % Se è la prima volta che incontriamo lo stato
                        % marcato in questo path, ci segnamo che lo
                        % abbiamo visitato e ci segnamo gli stati
                        % antecedenti ad esso
                        if (isequal(possible_next_states(i), marked_states(mar)) && ~visited_marked)
                            visited_marked = 1;
                            coac_states = union(coac_states, path);

                            % Altrimenti se abbiamo già visitato lo stato
                            % marcato e il prossimo stato di riferimento è
                            % di nuovo lo stato marcato ci segnamo il
                            % percorso attuale e aggiungiamo gli stati
                            % nuovi del path agli stati antecedenti al
                            % marcato (marcato non incluso)
                        elseif (isequal(possible_next_states(i), marked_states(mar)) && visited_marked)
                            all_paths{end+1} = path;
                            coac_states = union( coac_states, path( (index_of_marked+1):end ) );
                        end

                    end

                end

                % Finiamo quì se non ci sono transizioni successive per
                % questo stato o se tutte le transizioni successive sono
                % già state visitate (incluso lo stato marcato se già
                % visitato)
                % In entrambi i casi dovremmo resettare il path all'ultima
                % ramificazione valida
            else

                % Se è presente lo stato marcato nelle transizioni
                % successive allora ci segnamo il percorso attuale e
                % dimentichiamo tutti gli stati visitati segnati dopo
                % aver oltrepassato lo stato marcato, così che altri
                % path potranno arrivarci
                if any(strcmp(marked_states(mar), possible_next_states))
                    all_paths{end+1} = path;
                    index_from = find(strcmp( path(index_of_marked+1), all_visited_states), 1);
                    all_visited_states(index_from:end) = [];
                end

                if last_branch_index > 0
                    check_branch_value = transition_table( path(last_branch_index) );
                    successors = check_branch_value("successors");

                    num_of_next = size(successors,1);
                    temp = strings(1, 0);
                    tt = 1;
                    for nn = 1:num_of_next
                        row = successors{nn,2};
                        for mm = 1:length(row)
                            temp(tt) = row{mm};
                            tt = tt+1;
                        end
                    end
                    possible_next_states = temp;

                    % Se non abbiamo visitato lo stato marcato in
                    % questo path e tutte le transizioni successive
                    % dall'ultima ramificazione sono già state visitate
                    % possiamo ignorarla e passare direttamente alla
                    % ramificazione antecedente
                    if all(ismember(possible_next_states, temp_visited)) && ~visited_marked
                        branches(end) = [];
                        last_branch_index = branches(end);
                    end
                end

                % Reset del path all'ultima ramificazione valida
                path = path(1:last_branch_index);

                % Se abbiamo visitato lo stato marcato e stiamo
                % resettando il path ad un punto antecedente ad esso
                % allora possiamo ricominciare a tracciare le visite
                % sullo stato marcato da 0
                if last_branch_index < index_of_marked && visited_marked
                    visited_marked = 0;
                end

            end

        end
    end
end


%%% COMBINAZIONE DEI PATH

% Andiamo ora a considerare tutte le combinazioni di unione di path
% ricavati
num_of_paths = length(all_paths);

num_of_bits_needed = floor(log2(num_of_paths)) + 1;
for ii = 1:num_of_paths
    %combination = ;
end

% Non consideriamo la riga costituita da tutti zeri, quindi...
num_of_path_combinations = 2^num_of_paths - 1;
all_path_combination_matrix = dec2bin(0:2^(num_of_paths)-1) - '0';
all_path_combination_matrix(1,:) = [];

all_states_combinations = {};

for i = 1:num_of_path_combinations

    combination = all_path_combination_matrix(i,:);
    % paths_combination contiene solo i path considerati dalla
    % combinazione attuale
    paths_combination = all_paths(logical(combination));
    new_path = strings(0);

    % Uniamo tutti i path considerati dalla combinazione
    for j = 1:length(paths_combination)
        new_path = union(new_path, paths_combination{j});
    end

    % La funzione union restituisce un vettore colonna quindi facciamo
    % la trasposta per convenienza (potrebbe essere non necessario, ma
    % almeno aiuta nel debugging)
    new_path = new_path';

    if isempty(all_states_combinations)
        all_states_combinations{1} = new_path;

    else

        % Controllo per verificare che il nuovo path non è identico ad
        % un altro inserito in precedenza
        % Se risulta diverso, lo aggiungiamo alle combinazioni totali
        equal_check = 0;
        for k = 1:length(all_states_combinations)
            other_path = all_states_combinations{k};
            if isequal(sort(other_path), sort(new_path))
                equal_check = 1;
                break;
            end
        end
        if ~equal_check
            all_states_combinations{end+1} = new_path;
        end

    end
end


%%% CREAZIONE DEI SISTEMI DI TRANSIZIONE BASATI SUI PATH COMBINATI

% Creazione dei transition system basati su quello originale con
% pulizia delle matrici degli stati successivi/predecenti
all_sub_transition_systems = {};
all_states = keys(transition_table);
in_index = find(strcmp("inputs", all_states));
all_states(in_index) = [];

for j = 1:length(all_states_combinations)

    current_states = all_states_combinations{j};
    % Stati non considerati dalla combinazione attuale
    non_relevant_states = setxor(all_states, current_states);
    new_transition_table = containers.Map;
    new_transition_table("inputs") = transition_table("inputs");

    for i = 1:length(current_states)

        to_move = transition_table(current_states(i));
        values = containers.Map;
        values("initial") = to_move("initial");
        values("marked") = to_move("marked");
        values("output") = to_move("output");
        succ = to_move("successors");

        size_succ = size(succ,1);
        for stm = 1:size_succ
            stm_b = size_succ - stm + 1;
            row = succ{stm_b,2};
            len_row = length(row);
            for tf = 1:len_row
                tf_b = len_row - tf + 1;
                row_str = row{tf_b};
                if any(strcmp(row_str, non_relevant_states))
                    row(tf_b) = [];
                end
            end
            if isempty(row)
                succ(stm_b,:) = [];
            else
                succ{stm_b,2} = row;
            end
        end

        values("successors") = succ;

        pred = to_move("predecessors");

        size_pred = size(pred,1);
        for stm = 1:size_pred
            stm_b = size_pred - stm + 1;
            row = pred{stm_b,2};
            len_row = length(row);
            for tf = 1:len_row
                tf_b = len_row - tf + 1;
                row_str = row{tf_b};
                if any(strcmp(row_str, non_relevant_states))
                    row(tf_b) = [];
                end
            end
            if isempty(row)
                pred(stm_b,:) = [];
            else
                pred{stm_b,2} = row;
            end
        end

        values("predecessors") = pred;
        new_transition_table(current_states(i)) = values;

    end

    all_sub_transition_systems{end+1} = new_transition_table;

end


%%% CREAZIONE DELLE VARIAZIONI DEI SISTEMI DI TRANSIZIONE DOVUTE AI
%%% SELF-LOOP DEGLI STATI

self_loop_transition_systems_variations = {};

for ast = 1:length(all_sub_transition_systems)

    current_sts = all_sub_transition_systems{ast};
    all_states = keys(current_sts);
    in_index = find(strcmp("inputs", all_states));
    all_states(in_index) = [];

    % Matrice che tiene conto di tutti gli stati che hanno self-loop e
    % dei loro indici delle proprie matrici dei successori/predecessori
    all_selfloop_matrix = [];
    a = 1;

    for s = 1:length(all_states)

        values = current_sts(all_states{s});
        successors = values("successors");
        predecessors = values("predecessors");
        num_of_succ = size(successors,1);
        num_of_pred = size(predecessors,1);

        for r = 1:num_of_succ

            row_s = successors{r,2};
            row_s_div = [row_s{:}];
            self_index_s = find(strcmp(row_s_div, all_states{s}), 1);

            % Se abbiamo un self-loop, lo aggiungiamo alla matrice
            % insieme ai dati per ritrovarlo
            if ~isempty(self_index_s)

                for pr = 1:num_of_pred

                    row_p = predecessors{pr, 2};
                    row_p_div = [row_p{:}];
                    self_index_p = find(strcmp(all_states{s}, row_p_div), 1);
                    if ~isempty(self_index_p)
                        all_selfloop_matrix(a,1:5) = [s, r, self_index_s, pr, self_index_p];
                        a = a+1;
                        break;
                    end
                end


            end

        end

    end

    % Numero totale di self-loop nel sistema di transizione
    num_of_loops = size(all_selfloop_matrix,1);

    if num_of_loops > 0

        % Matrice di combinazione della presenza/assenza dei self-loop
        selfloop_combination_matrix = dec2bin(0:2^(num_of_loops)-1) - '0';

        for c = 1:size(selfloop_combination_matrix,1)

            % Combinazione che indica i self-loop da rimuovere
            combination = selfloop_combination_matrix(c,:);

            % Rimuoviamo i self-loop non considerati dalla combinazione
            % per la rimozione
            cut_selfloop_matrix = all_selfloop_matrix;
            cut_selfloop_matrix(~logical(combination), :) = [];

            % Creazione e ricostruzione di un sistema di transizione
            % basato sulla combinazione di self-loop attuale
            new_transition_system = containers.Map;

            for l = 1:length(all_states)

                og_value = current_sts(all_states{l});
                og_output = og_value("output");
                og_initial = og_value("initial");
                og_marked = og_value("marked");
                og_successors = og_value("successors");
                og_predecessors = og_value("predecessors");

                % Se lo stato corrente è presente nella matrice
                % self-loop
                if ~isempty(find(cut_selfloop_matrix(:,1)==l, 1))
                    for cc = 1:size(cut_selfloop_matrix,1)

                        % Rimozione dei self-loop di interesse
                        if cut_selfloop_matrix(cc,1) == l
                            og_succ_row = og_successors{all_selfloop_matrix(cc, 2),2};
                            og_succ_row(all_selfloop_matrix(cc, 3)) = [];
                            if isempty(og_succ_row)
                                og_successors(all_selfloop_matrix(cc, 2),:) = [];
                            else
                                og_successors{all_selfloop_matrix(cc, 2),2} = og_succ_row;
                            end
                            og_pred_row = og_predecessors{all_selfloop_matrix(cc, 4),2};
                            og_pred_row(all_selfloop_matrix(cc, 5)) = [];
                            if isempty(og_pred_row)
                                og_predecessors(all_selfloop_matrix(cc, 2),:) = [];
                            else
                                og_predecessors{all_selfloop_matrix(cc, 2),2} = og_pred_row;
                            end
                            %og_successors(all_selfloop_matrix(cc, 2), all_selfloop_matrix(cc, 3)) = [];
                            %og_predecessors(all_selfloop_matrix(cc, 4), all_selfloop_matrix(cc, 5)) = [];
                        end
                    end
                end

                new_value = containers.Map;
                new_value("initial") = og_initial;
                new_value("marked") = og_marked;
                new_value("output") = og_output;
                new_value("successors") = og_successors;
                new_value("predecessors") = og_predecessors;
                new_transition_system(all_states{l}) = new_value;

            end

            self_loop_transition_systems_variations{end+1} = new_transition_system;

        end

    else
        self_loop_transition_systems_variations{end+1} = current_sts;
    end

end

%%% RESTITUZIONE DELLA COLLEZIONE DI SOTTO-SISTEMI DI TRANSIZIONE

result = self_loop_transition_systems_variations;

end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Funzione per ottenere la parte accessibile di un transition system
function result = Ac(transition_table, initial_states)

%tic;
accessible_states = strings(0);
s = 1;
% Performiamo una ricerca DTS sul sistema di transizione a partire
% dagli stati iniziali
% Creamo uno stack dove inseriremo gli stati da visitare, andando
% dall'ultimo elemento al primo, aggiungendo stati successivi nuovi ed
% eliminando quelli già attraversati
stack = initial_states;
finished = 0;
while ~finished
    if isempty(stack) || ~isKey(transition_table, stack(end))
        finished = 1;
    else
        value = transition_table( stack(end) );
        possible_next_states = value("successors");

        num_of_next = size(possible_next_states,1);
        temp = strings(1, 0);
        tt = 1;
        for nn = 1:num_of_next
            row = possible_next_states{nn,2};
            for mm = 1:length(row)
                temp(tt) = row{mm};
                tt = tt+1;
            end
        end
        possible_next_states = temp;

        if ~any(strcmp(stack(end), accessible_states))
            accessible_states(s) = stack(end);
            s = s+1;
        end
        stack(end) = [];
        for l = 1:length(possible_next_states)

            state = possible_next_states(l);

            if ~any(strcmp(state, accessible_states))
                stack(end+1) = state;
            end

        end

    end

end

all_states = keys(transition_table);
in_index = find(strcmp("inputs", all_states));
all_states(in_index) = [];
non_accessible_states = setxor(all_states, accessible_states);
new_transition_table = containers.Map;
new_transition_table("inputs") = transition_table("inputs");
for i = 1:length(all_states)
    if any(strcmp(all_states{i}, accessible_states))
        %remove(transition_table, all_states{i});
        to_move = transition_table(all_states{i});
        values = containers.Map;
        values("initial") = to_move("initial");
        values("marked") = to_move("marked");
        values("output") = to_move("output");
        succ = to_move("successors");
        size_succ = size(succ,1);
        for stm = 1:size_succ
            stm_b = size_succ - stm + 1;
            row = succ{stm_b,2};
            len_row = length(row);
            for tf = 1:len_row
                tf_b = len_row - tf + 1;
                row_str = row{tf_b};
                if any(strcmp(row_str, non_accessible_states))
                    row(tf_b) = [];
                end
            end
            if isempty(row)
                succ(stm_b,:) = [];
            else
                succ{stm_b,2} = row;
            end
        end
        values("successors") = succ;

        pred = to_move("predecessors");
        size_pred = size(pred,1);
        for stm = 1:size_pred
            stm_b = size_pred - stm + 1;
            row = pred{stm_b,2};
            len_row = length(row);
            for tf = 1:len_row
                tf_b = len_row - tf + 1;
                row_str = row{tf_b};
                if any(strcmp(row_str, non_accessible_states))
                    row(tf_b) = [];
                end
            end
            if isempty(row)
                pred(stm_b,:) = [];
            else
                pred{stm_b,2} = row;
            end
        end
        values("predecessors") = pred;

        new_transition_table(all_states{i}) = values;
    end
end

result = new_transition_table;
%toc;

end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Funzione per ottenere la parte co-accessibile di un transition system
function result = CoAc(transition_table, marked_states)

%tic;
coaccessible_states = strings(0);
s = 1;
% Performiamo una ricerca DTS sul sistema di transizione a partire
% dagli stati marcati
% Creamo uno stack dove inseriremo gli stati da visitare, andando
% dall'ultimo elemento al primo, aggiungendo stati successivi nuovi ed
% eliminando quelli già attraversati
stack = marked_states;
finished = 0;
while ~finished
    if isempty(stack) || ~isKey(transition_table, stack(end))
        finished = 1;
    else
        value = transition_table( stack(end) );
        possible_previous_states = value("predecessors");

        num_of_prev = size(possible_previous_states,1);
        temp = strings(1, 0);
        tt = 1;
        for nn = 1:num_of_prev
            row = possible_previous_states{nn,2};
            for mm = 1:length(row)
                temp(tt) = row{mm};
                tt = tt+1;
            end
        end
        possible_previous_states = temp;

        if ~any(strcmp(stack(end), coaccessible_states))
            coaccessible_states(s) = stack(end);
            s = s+1;
        end
        stack(end) = [];
        for l = 1:length(possible_previous_states)

            state = possible_previous_states(l);

            if ~any(strcmp(state, coaccessible_states))
                stack(end+1) = state;
            end

        end

    end

end

all_states = keys(transition_table);
in_index = find(strcmp("inputs", all_states));
all_states(in_index) = [];
non_coaccessible_states = setxor(all_states, coaccessible_states);
new_transition_table = containers.Map;
new_transition_table("inputs") = transition_table("inputs");
for i = 1:length(all_states)
    if any(strcmp(all_states{i}, coaccessible_states))
        to_move = transition_table(all_states{i});
        values = containers.Map;
        values("initial") = to_move("initial");
        values("marked") = to_move("marked");
        values("output") = to_move("output");
        succ = to_move("successors");
        size_succ = size(succ,1);
        for stm = 1:size_succ
            stm_b = size_succ - stm + 1;
            row = succ{stm_b,2};
            len_row = length(row);
            for tf = 1:len_row
                tf_b = len_row - tf + 1;
                row_str = row{tf_b};
                if any(strcmp(row_str, non_coaccessible_states))
                    row(tf_b) = [];
                end
            end
            if isempty(row)
                succ(stm_b,:) = [];
            else
                succ{stm_b,2} = row;
            end
        end
        values("successors") = succ;

        pred = to_move("predecessors");
        size_pred = size(pred,1);
        for stm = 1:size_pred
            stm_b = size_pred - stm + 1;
            row = pred{stm_b,2};
            len_row = length(row);
            for tf = 1:len_row
                tf_b = len_row - tf + 1;
                row_str = row{tf_b};
                if any(strcmp(row_str, non_coaccessible_states))
                    row(tf_b) = [];
                end
            end
            if isempty(row)
                pred(stm_b,:) = [];
            else
                pred{stm_b,2} = row;
            end
        end
        values("predecessors") = pred;

        new_transition_table(all_states{i}) = values;
    end
end

result = new_transition_table;
%toc;

end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Funzione per ottenere la parte di un transition system contenente solo le
% traiettorie che vanno dagli stati iniziali agli stati marcati
function result = Trim(transition_table, initial_states, marked_states)

result = CoAc(Ac(transition_table, initial_states), marked_states);

end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~