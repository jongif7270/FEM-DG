f = @(x) sin(2*x)+cos(x);

t=-1:0.1:1;
y=f(t);

int_t=-1:0.5:1;
int_y=f(int_t);

num_int=length(int_t);

van=fliplr(vander(int_t));

w=van\int_y';

poly_int=@(x) x.^(0:num_int-1)*w;

z=zeros(1,length(t));

for i=1:length(t)
    z(i)=poly_int(t(i));
end


int_z=zeros(1,length(int_t));

for i=1:length(int_t)
    int_z(i)=poly_int(int_t(i));
end

figure
plot(t,y)
hold on
plot(int_t,int_y,'o')
hold on
%plot(int_t,int_z)
hold on
plot(t,z,'x')