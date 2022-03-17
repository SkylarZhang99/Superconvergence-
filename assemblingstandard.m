
function A = assemblingstandard(node,elem)
N=size(node,1); NT=size(elem,1);
A=zeros(N,N); %A = sparse(N,N);
for t=1:NT
     At=locatstiffness(node(elem(t,:),:));
     for i=1:3
         for j=1:3
        A(elem(t,i),elem(t,j))=A(elem(t,i),elem(t,j))+At(i,j);
         end
     end
end
