%--------------------------------------------------
%Note: This function solves the tridiagonal
%      system of equations using the thomas
%      algorithm.
%--------------------------------------------------
%http://books.google.com/books?id=7vuNLcQhg8UC&pg=PA42&lpg=PA42&dq=numerical+recipes+in+fortran+77+tridiagonal+solves&source=bl&ots=BCYH2EGCQH&sig=4XzlY3B7zImgBAux_7om1P8J3-8&hl=en&sa=X&ei=TX1SVLL_A7OasQSv7YD4Cg&ved=0CCIQ6AEwAA#v=onepage&q=numerical%20recipes%20in%20fortran%2077%20tridiagonal%20solves&f=false
function [q] = thomas(aa,bb,cc,dd)

bet(1)=bb(1);
gam(1)=dd(1)/bb(1); // u(1) = r(1)/bet
n=length(bb);

for i=2:n
    bet(i)=bb(i)-(aa(i)*cc(i-1)/bet(i-1)); // gam(i) = c(j-1)/bet(1)
    gam(i)=(dd(i)-aa(i)*gam(i-1))/bet(i);
end 

q(n)=gam(n);

for i=n-1:-1:1
   q(i)=gam(i)-(cc(i)*q(i+1)/bet(i));   
end 
