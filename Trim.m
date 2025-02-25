% Funzione per ottenere la parte di un transition system contenente solo le
% traiettorie che vanno dagli stati iniziali agli stati marcati
function result = Trim(transition_table, initial_states, marked_states)

    result = CoAc(Ac(transition_table, initial_states), marked_states);

end