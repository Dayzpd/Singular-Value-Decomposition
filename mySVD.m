%Zachary Day/Curtis Wesley
%December 26, 2017
%This functions takes in the mxn matrix A and the number of iterations for
%the QR algorithm

function [U,S,V] = mySVD(A,iter)
    %Creates a nxn matrix P using the matrix A. P is used to find the 
    %singular values of A. 
    P = A' * A;

    % We need to know the size of A to correctly compute U, V and S. 
    m=size(A,1);
    n=size(A,2);

    %Define another copy of the matrix P, will be used later to calculate
    %columns of V. 
    H=P;

    %QR algorithm to determine the eigenvalues of P. The eigenvalues are along
    %the main diagonal of the matrix. 
    for i=1:iter
        [Q,R] = qr(P); 
        P=R*Q; 
    end 

    myeig = [];
    myeig=double(myeig);
    for i=1:n
        myeig(i)=P(i,i);
    end
    
    %Create the nxn matrix V. It will be composed of the eigenvectors of the
    %matrix P. This can be found by computing the nullspace of (A-lambda*I).
    V = [];
    V=double(V);
    %Here is where we use the unmodified copy of P, denoted H. 
    for i=1:n
        AlamI = H - myeig(i) * eye(n);
        nullspace=null(AlamI); 
        V=[V (nullspace/norm(nullspace))];
    end
    %This will quit the program if QR method has not found all of the
    %eigenvalues. Convergence can be slow. Variance factors may determine
    %rate of convergence and it is heavily dependent on the matrix itself.
    
    if size(V,2)<size(V,1)
        disp('QR algorithm failed, please increase the number of iterations of algorithm.');
        quit cancel;
    end 
    
    %Creates a mxn generalized diagonal matrix S. The matrix S is not
    %necessarily square or invertible and therefore does not always have an inverse!
    S=zeros(m,n);
    minS=min(m,n);
    for i=1:minS
        S(i,i)=sqrt(myeig(i));
    end
    
  
    %Creates a mxm matrix U. To compute the columns of U, use the algorithm
    %Ui=A*Vi/sii, where Ui is the ith column of U, Vi is the ith column of
    %V and sii is the ith singular value and A is the original matrix. If
    %the number of columns of U is greater than the singular values, then
    %the rest of the columns of U are zero except there are 1's along the
%     %main diagonal. 
    U=eye(m);
    for i=1:n
        U(:,i) = A * V(:,i)/S(i,i);
    end
    size(A)
    size(U)
    size(S)
    size(V)
 
end