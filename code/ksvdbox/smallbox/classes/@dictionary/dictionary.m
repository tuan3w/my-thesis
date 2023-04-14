classdef dictionary
    %% Dictionary for sparse representation
    properties
        phi         %Matrix containing the dictionary
        len         %Length of basis functions
        nAtoms      %Number of basis function
        name        %String containing the matrix ensemble from which the dictionary is drawn
    end
    properties (Dependent = true)
        redundancy      %Redundancy of the dictionary: nAtoms/len
        coherence       %Maximum inner product of different basis
        isNormalised    %True if the atoms have unit norm
        rank            %rank of the dictionary
    end
    
    methods
        %% Constructor
        function obj = dictionary(phi,len,nAtoms)
            % obj = dictionary(phi,len,nAtoms)
            % INPUTS:
            % - phi: either a string specifying a matrix ensamble or a
            % matrix defining an explicit dictionary
            % - len: length of the atoms (only for implicit dictionaries)
            % - nAtoms: number of atoms (only for implicit dictionaries)
            if nargin
                if ~ischar(phi)
                    [obj.len obj.nAtoms] = size(phi);
                    obj.phi              = phi;
                    obj.name             = 'explicit';
                else
                    switch lower(phi)
                        case 'dct'
                            obj.phi = dctmatrix(len,nAtoms);
                        case 'grassmannian'
                            obj.phi = grassmanian(len,nAtoms);
                        otherwise
                            obj.phi = MatrixEnsemble(len,nAtoms,phi);
                    end
                    obj.len    = len;
                    obj.nAtoms = nAtoms;
                    obj.name   = lower(phi);
                end
            end
		end
        %% Dependent properties
        function redundancy = get.redundancy(obj)
            redundancy = obj.nAtoms/obj.len;
        end
        function coherence = get.coherence(obj)
            obj.phi = normcols(obj.phi);
            G = obj.phi'*obj.phi;
            G = G - eye(size(G));
            coherence = max(abs(G(:)));
        end
        function isNormalised = get.isNormalised(obj)
            isNormalised = norm(sum(conj(obj.phi).*obj.phi) - ...
                ones(1,obj.nAtoms))<1e-9;
        end
        function r = get.rank(obj)
            r = rank(obj.phi);
        end
        %% Operations
        function obj = normalize(obj)
            obj.phi = normcols(obj.phi);
		end
        %% Visualization
        function image(obj)
            %Image of the dictionary
            if isreal(obj.phi)
                imagesc(obj.phi);
                title('Dictionary');
                xlabel('Atom number');
            else
                subplot(2,1,1)
                imagesc(real(obj.phi));
                title('Real');
                xlabel('Atom number');
                subplot(2,1,2)
                imagesc(imag(obj.phi));
                title('Imaginary');
                xlabel('Atom number');
            end
        end
        function imagegram(obj)
            G = obj.phi'*obj.phi;
            imagesc(G);
            title('Gram Matrix')
        end
        function plot(obj,n)
            %Plot of the n-th basis
            if isreal(obj.phi)
                plot(obj.phi(:,n));
                title(['Atom number ' num2str(n) '/' num2str(size(obj.phi,2))]);
            else
                subplot(2,1,1)
                plot(real(obj.phi(:,n)));
                title(['Atom number ' num2str(n) '/' num2str(size(obj.phi,2)) ' - Real']);
                subplot(2,1,2)
                plot(imag(obj.phi(:,n)));
                title(['Atom number ' num2str(n) '/' num2str(size(obj.phi,2)) ' - Imaginary']);
            end
		end
        function movie(obj)
            %Movie of the basis
            for i=1:size(obj.phi,2)
                obj.plot(i);
                pause(1/25);
            end
        end
    end
end
