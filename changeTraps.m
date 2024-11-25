  %% draw random traps position
  function traps=changeTraps(trap_dist,len)
    trapvector_x=linspace(-trap_dist.*(floor(len/2)),trap_dist.*(floor(len/2)),len);
    trapvector_y=linspace(-trap_dist.*(floor(len/2)),trap_dist.*(floor(len/2)),len);
    trapmatrix_x=repmat(trapvector_x,len,1);
    trapmatrix_y=(repmat(trapvector_y,len,1))';
    move_x=(2.*trap_dist).*rand(len,len)-trap_dist;
    move_y=(2.*trap_dist).*rand(len,len)-trap_dist;
    trapmatrix_x=trapmatrix_x+move_x;
    trapmatrix_y=trapmatrix_y+move_y;
    traps(:,1)=trapmatrix_x(:);
    traps(:,2)=trapmatrix_y(:);
    cfg.initTrapPositions = traps;
