function D = dctmatrix(N,K,type)

error(nargchk(1,3,nargin,'struct'));
if ~exist('type','var') || isempty(type), type='II'; end
if ~exist('K','var') || isempty(K), K=N; end

[c r] = meshgrid(0:K-1,0:N-1);
switch type
    case 'I'
        D = cos(pi*c.*r/(K-1));
        D(1,:) = D(1,:)/sqrt(2);
        D(end,:) = D(end,:)/sqrt(2);
    case 'II'
        D = cos(pi*(2*c+1).*r/(2*K));
        D(1,:) = D(1,:)/sqrt(2);
    case 'III'
        D = cos(pi*(2*r+1).*c/(2*K));
        D(:,1) = D(:,1)/sqrt(2);
    case 'IV'
        D = cos(pi*(2*r+1+2*c+4*c.*r)/(4*K));
    otherwise
        error('unsupported dct type');
end
D = normcols(D);
