function [D,err,err2] = nsdl(params,varargin)

%%%%% parse input parameters %%%%%


data = params.data;

ompparams = {'checkdict','off'};

% iteration count %

if (isfield(params,'iternum'))
  iternum = params.iternum;
else
  iternum = 10;
end

% status messages %

printiter = 0;
printreplaced = 0;
printerr = 1;
printgerr = 0;

verbose = 't';
msgdelta = -1;

for i = 1:length(varargin)
  if (ischar(varargin{i}))
    verbose = varargin{i};
  elseif (isnumeric(varargin{i}))
    msgdelta = varargin{i};
  else
    error('Invalid call syntax');
  end
end

for i = 1:length(verbose)
  switch lower(verbose(i))
    case 'i'
      printiter = 1;
    case 'r'
      printiter = 1;
      printreplaced = 1;
    case 't'
      printiter = 1;
      printerr = 1;
      if (isfield(params,'testdata'))
        printgerr = 1;
      end
  end
end

if (msgdelta<=0 || isempty(verbose))
  msgdelta = -1;
end

ompparams{end+1} = 'messages';
ompparams{end+1} = msgdelta;



% compute error flag %

comperr = (nargout>=3 || printerr);


% validation flag %

testgen = 0;
if (isfield(params,'testdata'))
  testdata = params.testdata;
  if (nargout>=4 || printgerr)
    testgen = 1;
  end
end


% data norms %



% mutual incoherence limit %

% determine dictionary size %

if (isfield(params,'dictsize'))    % this superceedes the size determined by initdict
  dictsize = params.dictsize;
end

%if (size(data,2) < dictsize)
  %error('Number of training signals is smaller than number of atoms to train');
%end

H = params.H;
%[l, ~] = find(H == 1);
idx = params.idx;
classdata = {};
nclass = size(H,1);
subdictsize = params.subdictsize;

for i = 1:nclass
    classdata{i} = data(:,idx{i});
end
% initialize the dictionary %


%printf('\nInit from data')
D= [];

for i=1:nclass
    %break
    data_ids = idx{i};
    data_i = data(:, data_ids);
      %dic_params.data= data_i;
      %D_i = ksvd(dic_params, '');
    %data_ids = find(colnorms_squared_new(data_i) > 1e-6);   % ensure no zero data elements are chosen
    perm = randperm(length(data_ids));
    D_i = data(:,data_ids(perm(1:subdictsize)));
    D = [D D_i];
end



% normalize the dictionary %

%D = rand(size(data,1), dictsize);
D = normcols(D);

err = zeros(1,iternum);
err2 = zeros(1,iternum);




%%%%%%%%%%%%%%%%%  main loop  %%%%%%%%%%%%%%%%%

%printf('\nPrecompute Y*Y^T');
lambda2 = params.lambda2;
offset = 1;

for iter = 1:iternum
    err(iter) = 0.0;
    err2(iter) = 0.0;
  
    offset=1;
    %G = [];
    p2 = randperm(nclass);
    for cl = 1:nclass
		
        %cl = p2(cl2);
        %printf('\nOtimize dict of class %d', cl)
        data_i = classdata{cl};
        dictidx = (cl -1 ) * subdictsize + 1 : cl * subdictsize;
        D_i = D(:, dictidx);
		%Gamma = lsqnonneg(D_i, data_i);
		Gamma = ((D_i'* D_i + lambda2 * eye(size(D_i,2))) \ D_i') * data_i;
        for id =1:subdictsize
			%break
            %id = p(j);
            dj = D_i(:,id);
			dj_old = dj;
            x_j = Gamma(id,:);
            dx = dj * x_j;
            E_i =  data_i - D_i * Gamma + dx;
            d_j = E_i * x_j';
            d_j = d_j./norm(d_j);
            D_i(:,id) = d_j;
            if (msgdelta>0)
                timereta(tid, j, msgdelta);
            end
        end
        if (comperr)
            Gamma = (D_i'* D_i + lambda2 * eye(size(D_i,2)))\ (D_i' * data_i);
            err(iter) = err(iter) + compute_err(D_i,Gamma,data_i) + lambda2 * sum(sum(Gamma.^2));
        end
        D(:, dictidx) = D_i;
        offset = offset + subdictsize;
    end
    if (msgdelta>0)
        printf('updating atoms: iteration %d/%d', dictsize, dictsize);
    end

    %%%%%  compute error  %%%%%

    if (testgen)
        if (memusage >= MEM_NORMAL)
            G = D'*D;
        end
        GammaG = sparsecode(testdata,D,XtXg,G,thresh);
        gerr(iter) = compute_err(D,GammaG,testdata);
    end


    %%%%%  print info  %%%%%

    info = sprintf('Iteration %d / %d complete', iter, iternum);
    if (printerr)
        info = sprintf('%s, %s = %.4g', info, 'RSME', err(iter));
    end
    if (printgerr)
        info = sprintf('%s, test %s = %.4g', info, 'RSME', gerr(iter));
    end
    if (printreplaced)
        info = sprintf('%s, replaced %d atoms', info, sum(replaced_atoms) + cleared_atoms);
    end

    if (printiter)
        disp(info);
        if (msgdelta>0), disp(' '); end
    end

end
end