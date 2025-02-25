% Funzione per ottenere la parte accessibile di un transition system
function result = Ac(transition_table, initial_states)

    
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
        if isempty(stack)
            finished = 1;
        elseif ~isKey(transition_table, stack(end))
            stack(end) = [];
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
    non_accessible_states = setxor(all_states, accessible_states);
    new_transition_table = containers.Map;
    
    for i = 1:length(all_states)
    
        if any(strcmp(all_states{i}, accessible_states))
    
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

end