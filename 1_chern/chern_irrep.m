%-----------------------------------------------------------------------------------
%   0_chern_irrep.ipynb
% 
%   Calculates the Chern number associated with a particular irrep
%
%   Code written by Canon Sun. Last updated 5th September, 2024.
%-----------------------------------------------------------------------------------

% The code takes in two sets of data:
% 1.) The Hamiltonian at various fluxes. 
load('HAll_12_L=10.mat')

% 2.) The projector
Ns = 16; %Number of sublattices
PReal = readmatrix('G24_projector_03_real.txt');
PImag = readmatrix('G24_projector_03_imag.txt');
P = PReal + 1j*PImag;
P = kron(P,eye(Ns));



%Get all eigenvectors that transform in irrep
[VSorted, EngIndex] = Eigenvectors(HAll,P); 

%Select out only filled states
f=11/16; %Filling: 5/16, 8/16, or 11/16 for gaps
VFilled = FilledBands(VSorted,EngIndex,f);

%Compute chern number
C = ChernNumber(VFilled)

function [VSorted,EngIndex] = Eigenvectors(HAll,P)
    %This function simultaneously diagonalizes the Hamiltonian and
    %projector at all flux points. I then selects out states that transform
    %according to the chosen irrep. THe output is `VSorted`, a cell array
    %containing the eigenvectors corresponding to the selected irrep at
    %each flux point, and `EngIndex`, an index that indivates the position of
    %these eigenvectors when all eigenvalues are sorted in ascending order.

    %Parameters
    L = size(HAll,1); %Size of flux grid
    N = size(HAll,3); %Size of Hamiltonian
    
    %Initialize output
    VSorted = cell(L,L);
    EngIndex = cell(L,L);
    
    for i = 1:L
        for j = 1:L
            %For each flux point, perform simultaneous diagonalization
            Ham = squeeze(HAll(i,j,:,:));
            [V,Eng,D] = simdiag(Ham,P);
            [~,ind] = sort(real(diag(Eng)));

            %Count the number of states in irrep
            counter = 0;
            for b = 1:N
                if abs(D(b,b))<1.001 && abs(D(b,b))>0.999
                   counter = counter+1;
                end
            end

            %Now collect all states in that irrep
            VIrrep = zeros(N,counter);
            EngIndexIrrep = zeros(counter,1);
               
            counterp = 1;
            for b = 1:N
                if abs(D(b,b))<1.001 && abs(D(b,b))>0.999
                   VIrrep(:,counterp) = V(:,b);
                   EngIndexIrrep(counterp) =  ind(b);
                   counterp = counterp+1;
                end
            end       

            %Save eigenvectors and indices to their corresponding cell
            VSorted{i,j} = VIrrep;
            EngIndex{i,j} = EngIndexIrrep;
        end
    end
end

function VFilled = FilledBands(VSorted,EngIndex,f)
    %Select out eigenvectors that are filled

    %Parameters
    L = size(VSorted,1); %Flux grid dimension
    N = size(VSorted{1,1},1); %Size of Hamiltonian
    Nf = N*f; %Number of filled states
    
    VFilled=cell(L,L);
    for i=1:L
        for j=1:L
            %Extract data from cell
            VIrrep=VSorted{i,j};
            EngIndexIrrep=EngIndex{i,j};
            NStates=size(EngIndexIrrep,1);

            %Count the number of filled states
            counter = 0;
            for b=1:NStates
                if EngIndexIrrep(b)<=Nf
                   counter = counter+1;
                end
            end

            %Now construct output
            VIrrepFilled = zeros(N,counter);
            counterp = 1;
            for b = 1:counter
                if EngIndexIrrep(b) <= Nf
                    VIrrepFilled(:,b) = VIrrep(:,b);
                    counterp = counterp+1;
                end
            end
            VFilled{i,j} = VIrrepFilled;
        end
    end   
end


function C=ChernNumber(V)
    %Computes Chern number
    L = size(V,1);%Dimension of flux space
  
    %Initialize output
    C = 0;
    %Linking variable
    U=zeros(L-1,L-1,2);
    for i=1:L-1
        for j=1:L-1
            inew = mod(i,L-1)+1;
            jnew = mod(j,L-1)+1;

            V1 = V{i,j};
            V2 = V{inew,j};
            V3 = V{i,jnew};
               
            U(i,j,1) = det(V1'*V2);
            U(i,j,1) = U(i,j,1)/abs(U(i,j,1));
            U(i,j,2) = det(V1'*V3);
            U(i,j,2) = U(i,j,2)/abs(U(i,j,2));                
        end
    end
    %Chern number
    for i = 1:L-1
       for j = 1:L-1
           inew = mod(i,L-1)+1;
           jnew = mod(j,L-1)+1;
           F(i,j) = log(U(i,j,1)*U(inew,j,2)*conj(U(i,jnew,1))*conj(U(i,j,2)));
           C = C+F(i,j);
       end
   end
   C = -C/(2*pi*1j);
end

