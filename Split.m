% Funzione per ottenere tutti i sotto-sistemi di transizione di un dato
% transition system
function result = Split(transition_table, initial_states, marked_states)

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
                    
                    % Se abbiamo già visitato lo stato marcato e
                    % stiamo andando su uno stato antecedente lo stato marcato
                    % già visitato ci segnamo il path
                    if ~any( strcmp(possible_next_states, marked_states(mar)) ) && visited_marked
                        all_paths{end+1} = path;
                    end
    
                    % Se è presente lo stato marcato nelle transizioni
                    % successive allora ci segnamo il percorso attuale e
                    % dimentichiamo tutti gli stati visitati segnati dopo
                    % aver oltrepassato lo stato marcato, così che altri
                    % path potranno arrivarci
                    if any(strcmp(marked_states(mar), possible_next_states)) && index_of_marked < length(path)
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
    
    for j = 1:length(all_states_combinations)
    
        current_states = all_states_combinations{j};
        % Stati non considerati dalla combinazione attuale
        non_relevant_states = setxor(all_states, current_states);
        new_transition_table = containers.Map;
    
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