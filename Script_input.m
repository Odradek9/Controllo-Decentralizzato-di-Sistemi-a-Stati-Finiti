close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                            INIZIALIZZAZIONE                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plants_states = cell(1,2);
plants_initial_states = cell(1,2);
plants_inputs = cell(1,2);
plants_state_transitions = cell(1,2);
plants_transition_tables = cell(1,2);
specification_states = cell(1,2);
specification_inputs = cell(1,2);
specification_initial_states = cell(1, 2);
specification_marked_states = cell(1, 2);
specification_transition_tables = cell(1, 2);
specification_transition_matrices = cell(1, 2);
specification_state_output_cell_matrices = cell(1, 2);
tolerances = zeros(1, 2);

for p = 1:2

    disp("----------------------------------------------------------------------")
    disp(" ")

    disp("<strong>For plant " + string(p) + " :</strong>");
    disp(" ");

    %----------------------------------------------------------------------

    % Chiede di inserire tutti gli stati del plant, Xi
    disp("Would you like to import the plant " + string(p) + " states from file (1) or input the data manually (2) ?")
    disp('(File import compatible with outputs of MATLAB function "writematrix" in .txt format)')
    disp('(One state (numerical values/arrays only) for each row of the matrix)')
    while true
        response = input("");

        disp(" ")

        if response == 1
            data = read_plant_states(input("File name: ", "s"));
            break
        elseif response == 2
            data = create_plant_states();            
            break
        else
            disp("Invalid input. Please provide either 1 or 2 as an answer.")
        end
    end

    pi_states = data;

    %%%
    plants_states{p} = pi_states;

    disp(" ");

    %----------------------------------------------------------------------

    % Chiede di inserire tutti gli stati inziali del plant, Xi0
    disp("Plant states:")
    disp(pi_states);

    pi_initial_states = zeros(0);

    disp("Please insert the initial plant states one at a time: ")
    disp("(Minimum state number is 1, must be numerical or numerical array. Empty input to stop)")
    i = 1;
    while true
        temp = input('');
        if isnumeric(temp)
            if isempty(temp) && length(pi_initial_states) >= 1
                break
            elseif isempty(temp) && length(pi_initial_states) < 1
                disp("(Minimum number is 1, empty input to stop)");
            else
                if length(temp) == size(pi_states,2) && ismember(temp, pi_states, "rows")
                    if isempty(pi_initial_states) || ~ismember(temp, pi_initial_states, "rows")
                        pi_initial_states(i,:) = temp;
                        i = i+1;
                    else
                        disp("(Input was already given)")
                    end
                else
                    disp("(State given is not a valid plant state. Please try again)")
                end
            end
        else
            disp("(Given data is not numeric. Please insert proper data.");
        end
    end

    %%%
    plants_initial_states{p} = pi_initial_states;

    disp(" ");

    %----------------------------------------------------------------------

    % Chiede di inserire tutti gli input del plant, Ui
    disp("Would you like to import the plant " + string(p) + " inputs from file (1) or input the data manually (2) ?")
    disp('(File import compatible with outputs of MATLAB function "writematrix" in .txt format)')
    disp('(One input (numerical values/arrays only) for each row of the matrix)')
    while true
        response = input("");

        disp(" ")

        if response == 1
            data = read_plant_inputs(input("File name: ", "s"));
            break
        elseif response == 2
            data = create_plant_inputs();            
            break
        else
            disp("Invalid input. Please provide either 1 or 2 as an answer.")
        end
    end

    pi_inputs = data;

    %%%
    plants_inputs{p} = pi_inputs;

    disp(" ");

    %----------------------------------------------------------------------

    disp("Would you like to import the plant " + string(p) + " state transition map from file (1) or input the data manually (2) ?")
    disp('(File import compatible with output of given function "flatten_Fi_cell")')
    disp('(File import also compatible with outputs of MATLAB function "writecell" in .txt format, applied to the "flattened" cell structure)')
    disp('(Each row of cell must have the following structure: {xi, x3-i, ui, xi(t+1)_1, xi(t+1)_2, ...} )')
    while true
        response = input("");

        disp(" ")

        if p == 1
            length_i = size(pi_states,2);
            while true
                temp = input("Please state the plant 2 state dimensions: ");
                if ~isnumeric(temp)
                    disp("(Value must be numerical)")
                else
                    break;
                end
            end
            length_3_i = temp;
        else
            length_i = size(pi_states,2);
            p3_i_states = plants_states{1};
            length_3_i = size(p3_i_states,2);
        end

        if response == 1
            data = read_Fi_cell(input("File name: ", "s"), length_i, length_3_i, size(pi_inputs,2));
            break
        elseif response == 2
            data = create_Fi_cell(pi_states, pi_inputs);
            break
        else
            disp("Invalid input. Please provide either 1 or 2 as an answer.")
        end
    end

    Fi_transition_cell_matrix = data;

    %%%
    plants_state_transitions{p} = Fi_transition_cell_matrix;

    disp(" ")

    %----------------------------------------------------------------------

    pi_transition_table = containers.Map;
    pi_transition_table("inputs") = pi_inputs;

    for ps = 1:size(pi_states,1)
        values = containers.Map;
        if ismember(pi_states(ps,:), pi_initial_states, "rows")
            values("initial") = 1;
        else
            values("initial") = 0;
        end
        succ_cell_matrix = {};
        i = 1;
        for kl = 1:size(Fi_transition_cell_matrix,1)
            transition_condition = Fi_transition_cell_matrix{kl,1};
            if isequal(transition_condition{1}, pi_states(ps,:))
                succ_cell_matrix{i,1} = Fi_transition_cell_matrix{kl,1};
                succ_cell_matrix{i,2} = Fi_transition_cell_matrix{kl,2};
                i = i+1;
            end
        end
        values("successors") = succ_cell_matrix;
        pi_transition_table(num2str(pi_states(ps,:))) = values;
    end

    %%%
    plants_transition_tables{p} = pi_transition_table;

    %----------------------------------------------------------------------
    %----------------------------------------------------------------------

    % Chiede l'inserimento degli stati Xqi in uno string array
    disp("Would you like to import the plant " + string(p) + " specification states from file (1) or input the data manually (2) ?")
    disp('(File import compatible with outputs of MATLAB function "writematrix" in .txt format)')
    disp('(String values only)')
    while true
        response = input("");

        disp(" ")

        if response == 1
            data = read_plant_spec_states(input("File name: ", "s"));
            break
        elseif response == 2
            data = create_plant_spec_states();
            break
        else
            disp("Invalid input. Please provide either 1 or 2 as an answer.")
        end
    end

    qi_states = data;

    %%%
    specification_states{p} = qi_states;

    disp(" ")

    %----------------------------------------------------------------------

    % Chiede di inserire tutti gli input della specifica, Uqi
    disp("Would you like to import the plant " + string(p) + " specification inputs from file (1) or input the data manually (2) ?")
    disp('(File import compatible with outputs of MATLAB function "writematrix" in .txt format)')
    disp('(String values only)')
    while true
        response = input("");

        disp(" ")

        if response == 1
            data = read_plant_spec_inputs(input("File name: ", "s"));
            break
        elseif response == 2
            data = create_plant_spec_inputs();
            break
        else
            disp("Invalid input. Please provide either 1 or 2 as an answer.")
        end
    end

    qi_inputs = data;
    
    %%%
    specification_inputs{p} = qi_inputs;

    %----------------------------------------------------------------------
    
    number_of_states_qi = length(qi_states);
    number_of_inputs_qi = length(qi_inputs);

    disp(" ");

    disp("Would you like to import the plant " + string(p) + " specification transition matrix from file (1) or input the data manually (2) ?")
    disp('(File import compatible with outputs of MATLAB function "writematrix" in .txt format)')
    while true
        response = input("");

        disp(" ")

        if response == 1
            data = read_transition_matrix_qi(input("File name: ", "s"), number_of_states_qi, number_of_inputs_qi);
            break
        elseif response == 2
            data = create_transition_matrix_qi(qi_states, qi_inputs);            
            break
        else
            disp("Invalid input. Please provide either 1 or 2 as an answer.")
        end
    end

    Mqi = data;

    %%%
    specification_transition_matrices{p} = Mqi;
    
    disp(" ");

    %----------------------------------------------------------------------
    
    % Chiede l'inserimento degli stati iniziali Xqi,0
    disp("Plant specification states:")
    disp(qi_states);

    qi_initial_states = strings(0); 
    i = 1;

    disp('Please insert the initial specification states one at a time:')
    disp("(Insert at least 1 initial state, blank insert to stop):");
    while true
        
        initial_state = input('', "s");

        if isempty(initial_state) && length(qi_initial_states) < 1
            disp("Please provide at least 1 valid state.")
        elseif isempty(initial_state) && length(qi_initial_states) >= 1
            break
        else
            if any(strcmp(qi_states, initial_state))
                qi_initial_states(i) = initial_state;
                i = i+1;
            else
                disp("Please provide a valid state.");
            end
        end
    end
    
    %%%
    specification_initial_states{p} = qi_initial_states;

    %----------------------------------------------------------------------
    
    % Chiede l'inserimento degli stati marcati Xi,q,m
    qi_marked_states = strings(0); 
    i = 1;
    disp('Please insert the marked specification states one at a time: ');
    disp("(Insert at least 1 marked state, blank insert to stop)");
    while true
        
        marked_state = input('', "s");
    
        if isempty(marked_state) && length(qi_marked_states) < 1
            disp("Please provide at least 1 valid state.")
        elseif isempty(marked_state) && length(qi_marked_states) >= 1
            break
        else
            if any(strcmp(qi_states, marked_state))
                qi_marked_states(i) = marked_state;
                i = i+1;
            else
                disp("Please provide a valid state.");
            end
        end
    end
    
    %%%
    specification_marked_states{p} = qi_marked_states;
    
    disp(" ");

    %----------------------------------------------------------------------
    
    % Chiede l'inserimento della relazione stato-output per ogni stato
    disp("Would you like to import the plant " + string(p) + " specification states outputs from file (1) or input the data manually (2) ?")
    disp('(File import compatible with outputs of MATLAB function "writematrix" in .txt format)')
    disp('(Numerical values/arrays only, must match previously given plant states)')
    disp("(Make sure the order of the outputs in the array matches the order of the plant specification states)")
    while true
        response = input("");

        disp(" ")

        if response == 1
            data = read_plant_spec_outputs(input("File name: ", "s"), pi_states);
            break
        elseif response == 2
            data = create_plant_spec_outputs(qi_states, pi_states);
            break
        else
            disp("Invalid input. Please provide either 1 or 2 as an answer.")
        end
    end

    output_vector = data;

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    % Crea la matrice stato-uscite utilizzando la funzione Hi,q
    Hqi = create_Hqi(qi_states, output_vector);
    %%%
    specification_state_output_cell_matrices{p} = Hqi;

    %----------------------------------------------------------------------
    
    % Creazione della mappa rappresentante il sistema di transizione
    
    specification_transition_table_Qi = containers.Map;
    specification_transition_table_Qi("inputs") = mat2cell(qi_inputs, 1, ones(1, numel(qi_inputs)));
    for i = 1:number_of_states_qi
        %complete_successor_vector = strings(number_of_inputs_qi,0);
        complete_successor_vector = cell(0,2);
        csv = 1;
        %complete_predecessor_vector = strings(number_of_inputs_qi,0);
        complete_predecessor_vector = cell(0,2);
        cpv = 1;
            
        for k = 1:number_of_inputs_qi
            partial_successor_vector = {};
            partial_predecessor_vector = {};
            suc = 1;
            pre = 1;
            
            for j = 1:number_of_states_qi

                if Mqi(i, j, k) == 1
                    %complete_successor_vector(k, suc) = qi_states(j);
                    partial_successor_vector{suc} = qi_states(j);
                    suc = suc + 1;
                end
                if Mqi(j, i, k) == 1
                    %complete_predecessor_vector(k, pre) = qi_states(j);
                    partial_predecessor_vector{pre} = qi_states(j);
                    pre = pre + 1;
                end
            end

            if ~isempty(partial_successor_vector)
                complete_successor_vector{csv,1} = qi_inputs(k);
                complete_successor_vector{csv,2} = partial_successor_vector;
                csv = csv + 1;
            end
            if ~isempty(partial_predecessor_vector)
                complete_predecessor_vector{cpv,1} = qi_inputs(k);
                complete_predecessor_vector{cpv,2} = partial_predecessor_vector;
                cpv = cpv + 1;
            end
            
        end

        marked = any(strcmp(qi_marked_states, qi_states(i)));
        initial = any(strcmp(qi_initial_states, qi_states(i)));
        output = Hqi{i ,2};
        
        value = containers.Map;
        value("initial") = initial;
        value("marked") = marked;
        value("output") = output;
        value("successors") = complete_successor_vector;
        value("predecessors") = complete_predecessor_vector;

        specification_transition_table_Qi(qi_states(i)) = value;
    end
    
    %%%
    specification_transition_tables{p} = specification_transition_table_Qi;
    
    disp(" ");

    %----------------------------------------------------------------------
    
    % Chiede di inserire il livello di tolleranza per le specifiche
    disp("Please insert the desired accuracy for the given specifications: ");
    tolerance_i = -1;
    while tolerance_i < 0
        temp = input("");
        if ~isnumeric(temp)
            disp("Invalid data, please insert a numeric value.");
        else
            tolerance_i = temp;
        end
    end
    
    %%%
    tolerances(p) = tolerance_i;

    %----------------------------------------------------------------------
    
    disp(" ");
    disp(" ");

end

clearvars -except plants_states plants_initial_states plants_inputs ...
    plants_state_transitions plants_transition_tables ...
    specification_initial_states specification_marked_states ...
    specification_transition_tables specification_transition_matrices ...
    specification_state_output_cell_matrices tolerances ...
    specification_inputs specification_states

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              FUNZIONI                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Funzione creazione matrice stato-uscita per transition system
function result = create_Hqi(q1_states, output_vector)
    
    state_output_cell_matrix = cell(length(q1_states), 2);
    for i = 1:length(q1_states)
        state_output_cell_matrix{i,1} = q1_states(i);
        state_output_cell_matrix{i,2} = output_vector(i,:);
    end

    result = state_output_cell_matrix;
    
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Funzione che prende da file .txt valori matrice automa qi
% Compatibile con format output di funzione writematrix()
function result = read_transition_matrix_qi(filename, number_of_states, number_of_inputs)

    while true
        matrix = readmatrix(filename);
        if isequal(number_of_states, size(matrix,1)) && isequal(number_of_inputs, size(matrix,2)/size(matrix,1))
            break
        else
            disp(" ")
            disp("Number of states and/or inputs previously given does not match matrix dimensions.");
            filename = input("Please provide another filename: ", "s");
        end
    end

    fixed_matrix = zeros(number_of_states, number_of_states, number_of_inputs);
    for i = 1:number_of_inputs
        from = (number_of_states*i - number_of_states+1);
        to = (i*number_of_states);
        fixed_matrix(:,:,i) = matrix( :, from:to );
    end
    
    disp("Data collected successfully.")
    
    result = fixed_matrix;

end

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

% Funzione che crea la matrice dell'automa delle specifiche da input utente
function result = create_transition_matrix_qi(states, inputs)
    
    disp(" ");

    number_of_states = length(states);

    number_of_inputs = length(inputs);

    % Inizializza matrice automa, definisce la relazione stato-stato-input 
    % per ogni stato
    automa_matrix = zeros(number_of_states, number_of_states, number_of_inputs);
    
    disp("All states: ");
    disp(states);
    
    % Chiede l'inserimento dei possibili stati successive dato un certo 
    % input, per ogni stato
    for i = 1:number_of_states
        for k = 1:number_of_inputs
            disp(" ");
            next_states = [];
            disp('For specification state "' + states(i) + '" , input "' + inputs(k) + '"');
            disp("please insert the possible transitions array :");
            disp("(Input NaN for empty transition)");
            while isempty(next_states)
                next_states = string(input(""));
        
                if ismatrix(next_states)
        
                    % Controllo se il vettore dat0 è 
                    % sotto-vettore del vettore degli stati
        
                    % Se abbiamo un sotto-vettore valido, andiamo a
                    % registrare la relazione di transizione nella matrice
                    % dell'automa
                    A = next_states;
    
                    % Controlla se una transizione dello stato i per input
                    % k esiste o meno
                    if ~ismissing(A)
                        B = intersect(A, states);
                        if size(A) == size(B)
                            for l = 1:length(A)
                                index = find(states==next_states(l));
                                automa_matrix(i, index, k) = 1;
                            end
                        else
                            disp("Invalid data, please make sure the states are valid.")
                            next_states = [];
                        end
                    end
            
                else
                    disp("Invalid data, please make sure the formatting is correct.")
                    next_states = [];
                end
            end
        end
    end

    result = automa_matrix;
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function result = create_plant_states()

    pi_states = zeros(0);
    disp("Please insert the plant states one at a time: ")
    disp("(Minimum state number is 2, must be numerical or numerical array. Empty input to stop)")
    i = 1;
    while true
        temp = input('');
        if isnumeric(temp)
            if isempty(temp) && length(pi_states) > 1
                break
            elseif isempty(temp) && length(pi_states) < 2
                disp("(Minimum number is 2, empty input to stop)");
            else
                if ~isempty(pi_states) && length(temp) ~= size(pi_states,2)
                    disp("(Given state does not align with previously given states)")
                elseif isempty(pi_states) || ~ismember(temp, pi_states, "rows")
                    pi_states(i,:) = temp;
                    i = i+1;
                else
                    disp("(State was already given)")
                end
            end
        else
            disp("(Given data is not numeric. Please insert proper data.");
        end
    end

    result = pi_states;

end

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

function result = read_plant_states(filename)

    matrix = readmatrix(filename);
    disp("Data collected successfully.")

    result = matrix;

end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function result = create_plant_inputs()

    % Chiede di inserire tutti gli input del plant, Ui
    pi_inputs = zeros(0);
    disp("Please insert the plant inputs one at a time: ")
    disp("(Minimum input number is 1, must be numerical or numerical array. Empty input to stop)")
    i = 1;
    while true
        temp = input('');
        if isnumeric(temp)
            if isempty(temp) && length(pi_inputs) >= 1
                break
            elseif isempty(temp) && length(pi_inputs) < 1
                disp("(Minimum number is 1, empty input to stop)");
            else
                if ~isempty(pi_inputs) && length(temp) ~= size(pi_inputs,2)
                    disp("(Given input does not align with previously given inputs)")
                elseif isempty(pi_inputs) || ~ismember(temp, pi_inputs, "rows")
                    pi_inputs(i,:) = temp;
                    i = i+1;
                else
                    disp("(Input was already given)")
                end
            end
        else
            disp("(Given data is not numeric. Please insert proper data.");
        end
    end

    result = pi_inputs;

end

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

function result = read_plant_inputs(filename)

    matrix = readmatrix(filename);
    disp("Data collected successfully.")

    result = matrix;

end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function result = flatten_Fi_cell(cell_to_flatten)

    flat_cell = cell(size(cell_to_flatten,1), 4);
    
    for i = 1:size(cell_to_flatten,1)
        inner_cell = cell_to_flatten{i,1};
        outputs = cell_to_flatten{i,2};

        flat_cell{i,1} = inner_cell{1};
        flat_cell{i,2} = inner_cell{2};
        flat_cell{i,3} = inner_cell{3};
        for j = 1:size(outputs,1)
            flat_cell{i,3+j} = outputs(j,:);
        end
        
    end

    result = flat_cell;
end

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

function result = unflatten_Fi_cell(cell_to_unflatten, state_i_length, state_3_i_length, input_length)

    original_cell = cell(size(cell_to_unflatten,1), 2);
    
    for i = 1:size(cell_to_unflatten,1)

        element_1 = zeros(1,state_i_length);
        in = 1;
        for k = 1:state_i_length
            element_1(1,in) = cell_to_unflatten{i,k};
            in = in+1;
        end
        inner_cell{1} = element_1;

        element_2 = zeros(1,state_i_length);
        in = 1;
        for k = (state_i_length+1):(state_i_length+state_3_i_length)
            element_2(1,in) = cell_to_unflatten{i,k};
            in = in+1;
        end
        inner_cell{2} = element_2;

        element_3 = zeros(1,input_length);
        in = 1;
        for k = (state_i_length+state_3_i_length+1):(state_i_length+state_3_i_length+input_length)
            element_3(1,in) = cell_to_unflatten{i,k};
            in = in+1;
        end
        inner_cell{3} = element_3;
        index_1 = state_i_length+state_3_i_length+input_length;
        outputs = [];
        for j = 1:( (size(cell_to_unflatten(i,:),2) - (state_i_length+state_3_i_length+input_length))/state_i_length )
            index_2 = index_1 + 1 + state_i_length*(j-1);
            if ismissing( cell_to_unflatten{i, 3 + 1 + state_i_length*(j-1) } )
                break;
            end
            element_4 = zeros(1,state_i_length);
            in = 1;
            for k = index_2:(index_2+state_i_length-1)
                element_4(1,in) = cell_to_unflatten{i, k};
                in = in+1;
            end
            outputs(j,:) = element_4;
        end
        original_cell{i,1} = inner_cell;
        original_cell{i,2} = outputs;
        
    end

    result = original_cell;
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function result = read_Fi_cell(filename, state_i_length, state_3_i_length, input_length)

    og_cell = readcell(filename);
    
    while true
        
        if size(og_cell, 2) >= 4
            break;
        else
            disp("File does not match the expected format. Cell array must be at least 4 cells wide.")
            filename = input("Please provide another filename: ", "s");
        end

    end

    disp("Data collected successfully.")

    result = unflatten_Fi_cell(og_cell, state_i_length, state_3_i_length, input_length);

end

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

function result = create_Fi_cell(pi_states, pi_inputs)

    Fi_transition_cell_matrix = cell(1,2);
    disp("Please input the plant F(xi,x3-i,ui) transitions:")
    disp("(Empty input if there is no transition)")
    disp(" ")
    i = 1;
    for ps = 1:size(pi_states, 1)
        disp("For plant state " + num2str(pi_states(ps, :)) + ", insert the transitions following the example:")
        disp("(<strong>Example: { x3-i, ui, [xi; xi; ...] } </strong>)")
        disp("(With x3-i being the other plant state, ui being an input and xi being the result vector.)")
        state_transition_counter = 0;
        while true
            temp = input('');
            if iscell(temp) || isempty(temp)
                if (isempty(temp)  && state_transition_counter > 0) || ~iscell(temp)
                    break;
                else
                    x3_i = temp{1};
                    ui = temp{2};
                    outputs = temp{3};
                    if ismember(ui, pi_inputs, "rows") && ( size(outputs, 1) == size(intersect(outputs, pi_states, "rows"), 1) )
                        Fi_transition_cell_matrix{i,1} = {pi_states(ps,:), x3_i, ui};
                        Fi_transition_cell_matrix{i,2} = outputs;
                        state_transition_counter = state_transition_counter +1;
                        i = i+1;
                        disp("- - - Transition accepted.")
                    else
                        disp("(One or more of the given data does not align with previously given data (plant states, inputs...). Please try again)")
                    end
                end
            else
                disp("(Given data is not a cell. Please insert proper data.");
            end
        end
    end

    result = Fi_transition_cell_matrix;

end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function result = read_plant_spec_states(filename)

    matrix = readmatrix(filename, 'OutputType', 'string');
    disp("Data collected successfully.")

    result = matrix;

end

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

function result = create_plant_spec_states()

    qi_states = strings(0);
    disp("Please insert the specification states, in order, one at a time:");
    disp("(Minimum number is 2, empty input to stop)");
    i = 1;
    while true
        disp('State <strong>' + string(i) + '</strong>: ');
        temp = string(input(''));
        if isempty(temp) && length(qi_states) > 1
            break
        elseif isempty(temp) && length(qi_states) < 2
            disp("(Minimum number is 2, empty input to stop)");
        elseif isequal(temp, "inputs")
            disp("(State with name 'inputs' cannot be assigned as it's a special key)")
        else
            qi_states(i) = temp;
            i = i+1;
        end
    end

    result = qi_states;

end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function result = read_plant_spec_inputs(filename)

    matrix = readmatrix(filename, 'OutputType', 'string');
    disp("Data collected successfully.")

    result = matrix;

end

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

function result = create_plant_spec_inputs()

    qi_inputs = strings(0);
    disp("Please insert the specification inputs one at a time: ")
    disp("(Minimum input number is 1. Empty input to stop)")
    i = 1;
    while true
        temp = string(input(''));
        
        if isempty(temp) && length(qi_inputs) >= 1
            break
        elseif isempty(temp) && length(qi_inputs) < 1
            disp("(Minimum number is 1, empty input to stop)");
        else
            if isempty(qi_inputs) || ~any(strcmp(temp, qi_inputs))
                qi_inputs(i) = temp;
                i = i+1;
            else
                disp("(Input was already given)")
            end
        end
        
    end

    result = qi_inputs;

end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function result = read_plant_spec_outputs(filename, pi_states)

    correct = false;
    
    while ~correct
        matrix = readmatrix(filename);

        correct = true;
        for i = 1:size(matrix, 1)
            if ~( ismember(matrix(i, :), pi_states, 'rows') || all(isnan(matrix(i,:))) )
                disp("(One or more output values do not match previously given plant states)");
                filename = input("Please provide another filename: ", "s");
                correct = false; 
                break;
            end
        end
    end
    
    disp("Data collected successfully.")

    result = matrix;

end

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

function result = create_plant_spec_outputs(qi_states, pi_states)

    number_of_states_qi = length(qi_states);

    output_vector = zeros(number_of_states_qi,size(pi_states,2));
    for i = 1:number_of_states_qi
        temp = NaN;
        while true
            temp = input("Please insert the output for state <strong>" + string(qi_states{i}) + "</strong> :");
            if isnumeric(temp)
                if ( length(temp) == size(pi_states,2) && ismember(temp, pi_states, "rows") ) || ismissing(temp)
                    break;
                else
                    disp("(Given output is not a valid plant state.)");
                end
            else
                disp("(Given data is not numeric. Please insert proper data.")
            end
        end
        output_vector(i,:) = temp;
    end

    result = output_vector;

end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~