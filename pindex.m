function id = pindex
  t = getCurrentTask();
  id = iff(isempty(t), 0, get(t,'ID'));
end
