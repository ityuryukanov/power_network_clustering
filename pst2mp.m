function mpw = pst2mp(bus, lin, gen, basmva, f0, acdc, debug)

if nargin<5
  f0 = 50;
end
if nargin<6
  acdc = 'ac';
end
if nargin<7
  debug = false;
end

nb = size(bus,1);
nl = size(lin,1);
bus(:,1) = round(bus(:,1));
gen(:,1:2) = round(gen(:,1:2));
lin(:,1:2) = round(lin(:,1:2));

if size(bus,2)<11
  bus(:,11) = 100*ones(nb,1);
end
if size(bus,2)<12
  bus(:,12) = -100*ones(nb,1);
end
if size(bus,2)<13
  bus(:,13) = ones(nb,1);
end
if size(bus,2)<15
  bus(:,14:15) = [-0.5*ones(nb,1), 1.5*ones(nb,1)];
end

[~, ord] = sort(bus(:,1), 'ascend');
bus = bus(ord,:);
[~, ord] = sort(gen(:,2), 'ascend');
gen = gen(ord,:);

bus_mp = [bus(:,1), bus(:,10), bus(:,6:7), bus(:,8:9), ones(nb,1), bus(:,2:3), bus(:,13), ones(nb,1), bus(:,14:15)];
bus_typ = bus(:,10);
bus_pq = find(bus_typ==3);
bus_sw = find(bus_typ==1);
bus_mp(bus_pq,2) = 1;
bus_mp(bus_sw,2) = 3;
bus_mp(:,3:6) = bus_mp(:,3:6)*basmva;

busSnom = accumarray(gen(:,2), gen(:,7));
[busGidx,~,busSnom] = find(busSnom);
ng = size(busGidx,1);
isgen = ismember(bus(:,1), gen(:,2));  %round-off errors are possible!
gen_mp = [busGidx, bus(isgen,4:5), bus(isgen,11:12), bus(isgen,2), busSnom,...
  ones(ng,1), 2*abs(bus(isgen,4)), -2*abs(bus(isgen,4))];
gen_mp(:,2:5) = gen_mp(:,2:5)*basmva;

% Include generation in PST at non-generator buses into MPW-loads
nongenbus = setdiff(bus(bus(:,4)~=0|bus(:,5)~=0,1), gen_mp(:,1));
[~,mpw_idx,~] = intersect(bus_mp(:,1), nongenbus,'stable');
[~,pst_idx,~] = intersect(bus(:,1), nongenbus,'stable');
bus_mp(mpw_idx,3:4) = bus_mp(mpw_idx,3:4) - bus(pst_idx,4:5)*basmva;

[DUMMY,pstbus,pstlin] = ...
  evalc('PST.LLoadflow(bus,lin,2e-13,52,0.5,1.5,1.0,''n'',2)');
flo_fr = 1:2:2*nl-1;
flo_to = 2:2:2*nl;
lin_id = pstlin(flo_fr,1);
lin_mp = [lin(lin_id,1:2),lin(lin_id,3:5),zeros(nl,3),lin(lin_id,6:7),...
   ones(nl,1),zeros(nl,2),pstlin(flo_fr,4:5),pstlin(flo_to,4:5)];

mpw = struct_mp(bus_mp,lin_mp,gen_mp,basmva);
mpw.f0 = f0;
mpOpt = mpoption('pf.tol',2e-13,'pf.nr.max_it',50,'pf.enforce_q_lims',0,...
  'out.all',0,'verbose',0,'out.suppress_detail',1);
if strcmp(acdc,'ac')
  mpw = runpf(mpw,mpOpt);
else
  mpw = rundcpf(mpw,mpOpt);
end  
if debug
  mpw.bus(:,8:9) = pstbus(:,2:3);  % initialize voltages
  mpw = runpf(mpw,mpOpt);
  assert( mean( abs(pstbus(:,2) - mpw.bus(:,8)) )<1e-4 );
  assert( mean( abs(pstbus(:,3) - mpw.bus(:,9)) )<5e-3 );
end

end