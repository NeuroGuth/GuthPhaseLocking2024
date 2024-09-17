function TG_save(fname, data)

% Saves variables from inside parfor loop.
% For that purpose this external function has to be called

save(fname, '-struct', 'data');

end