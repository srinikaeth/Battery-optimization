function [ dcopf_soln ] = dcopf_storage_network_v19( Y,T,p_buy,p_sell,L_buses,slack,thetaL,thetaU, ...
                          n_storage,pen,n_buses,n_homes, storage_bus,X,R,Q_mult,TP)

n = length(Y);
B = imag(Y);

n_rem = n_homes - n_storage;

dcopf_soln.u = zeros(n_storage,T);
dcopf_soln.b = zeros(n_storage,T);

cvx_begin
    variable theta(n_buses,T);
    variable PG(T);        % PG is only at slack bus, so only one variable for each hr
        
    variable u(pen,T);
    variable b(pen,T+1);
    
    expression g;
    expression g_buy;
    expression g_sell;
    expression buy_opt;
    expression sell_opt;
    expression u_bus;   
    expression buy_opt_1;
    expression sell_opt_1;
    expression V;
    expression Q;
    
    %%% Defining expressions for the objective function
%         g_buy = sum(L_opf + [u; zeros(n_rem,T)]); 
%         buy_1 = max(g_buy,0);
%         
%         g_sell = max(L_opf + [u; zeros(n_rem,T)],0);
%         p_sell2 = -p_sell;
%         
%         buy_opt = (p_buy)'*(buy_1');  
%         sell_opt = (p_sell2)'*(g_sell');
%         
%         buy_opt_1 = sum(buy_opt);
%         sell_opt_1 = sum(sell_opt);
        
    minimize((p_buy - p_sell)'*max(PG,0) + p_sell'*PG) 
    subject to
    
        %%% Storage Constraints   
                
        for i = 1:pen 
            %%% Initial condition for the battery charging state
            b(i,1) == 0;
            
            %%% Battery relation
            for v = 1:T
                b(i,v+1) == b(i,v) + u(i,v);
            end
            
            %%% Charging state limits
            zeros(1,T+1) <= b(i,:) <= 0.243*ones(1,T+1);
            
            %%% Charging rate limits
            -0.090*ones(1,T) <= u(i,:) <= 0.090*ones(1,T);
        end
        
        %%%%% Calculating u_bus for Power demand (PD)
       for i = 1:T
            u_bus(1:pen,i) = u(1:pen,i);
            u_bus(pen+1:n_buses,i) = zeros(n_buses-pen,1);
        end

        %%%%%% Network Constraints
                       
        %%%% Bus angle constraints        
         thetaL <= theta <= thetaU;
         
        %%%% DC Power Flow Equation
        for i = 1:T
            for v = 2:n_buses
               theta(v,i)*sum((B(v,:))) - B(v,:)*theta(:,i) == -L_buses(v,i) - u_bus(v,i) ;
            end
            theta(1,i)*sum((B(1,:))) - B(1,:)*theta(:,i) == PG(i) -L_buses(1,i) - u_bus(1,i);
            
            %%%%% Transformer limits
             -TP <= L_buses(2:123,i) + u_bus(2:123,i) <= TP;
             
             %%%%% Power balance
             PG(i) == sum(L_buses(:,i) + u_bus(:,i));
        end
        
        
           
        %%% Slack angle constraint
        theta(slack,:) == 0;
       
        %%%% Voltage constraint
        Q = Q_mult.*(-L_buses - u_bus);
        V = R*(-L_buses - u_bus) + X*Q + 1;
        
        V(slack,:) == 1;
        
        0.95 <= V <= 1.1;
        
cvx_end

dcopf_soln.u = u;    
dcopf_soln.b = b;
dcopf_soln.Cost = cvx_optval;   
dcopf_soln.theta = theta;
dcopf_soln.u_bus = u_bus;
dcopf_soln.PG = PG;
dcopf_soln.buy = buy_opt;
dcopf_soln.sell = sell_opt;

dcopf_soln.PF_opt = zeros(n);
for i=1:n
    for j=i+1:1:n
        dcopf_soln.PF_opt(i,j)=B(i,j)*(theta(i)-theta(j));
    end
end

end
