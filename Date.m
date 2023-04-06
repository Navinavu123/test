function varargout = main_gui(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @main_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end


function main_gui_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;

axes(handles.axes1)
cla;

guidata(hObject, handles);





function varargout = main_gui_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)

function edit1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)

function edit2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)

function edit3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)

function edit4_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)

function edit5_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pushbutton1_Callback(hObject, eventdata, handles)


global A N BSx BSy x y
% Define area 
A = str2double(get(handles.edit2,'String'));

% Define number of nodes
N = str2double(get(handles.edit1,'String'));

% Basestation cordinates
BSx = str2double(get(handles.edit3,'String'));
BSy = str2double(get(handles.edit4,'String'));


%%
axes(handles.axes1)
cla;
% Generate random coordinates
x = A*rand(1,N);
y = A*rand(1,N);

% Plot the network
plot(x,y,'k.');
hold on

% Plot basestation
plot(BSx,BSy,'rs','markerfacecolor','g','markeredgecolor','k','markersize',12)
text(BSx-3,BSy-5,'BS')


axis([0 max(A,BSx) 0 max(A,BSy)]);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)


global A N BSx BSy x y
% Initial energy
Einit = str2double(get(handles.edit5,'String')); % Joules

% Packat length
PL = 400;

% Probability of becoming CH
p = str2double(get(handles.edit6,'String'));

if get(handles.radiobutton1,'Value')==1

    
    msk = 971; % Any prime number

    % Generate memory for sensor nodes
    memory = struct;

    % --------------------------------------------
    % Plot the network
    plot(x,y,'k.');
    hold on

    % Plot basestation
    plot(BSx,BSy,'rs','markerfacecolor','g','markeredgecolor','k','markersize',12)
    text(BSx-3,BSy-5,'BS')


    axis([0 A+30 0 A+30]);

    % Setup
    % Transmitt data to all nodes
    for ii = 1:N

        % -------------------------------------
        % Store in sensor node memory
        memory(ii).masterkey = msk;
        % Show communication line
        plot([x(ii) BSx],[y(ii) BSy],'r-');

    end
    title('Setup')
    pause(0.1)

    % Generate prime numbers
    Ps = primes(1000);
    Ps(1:2) = [];
   
    for ii = 1:N % For all nodes

        % Get master key
        msk = memory(ii).masterkey;

        % Find out nearest prime numbers
        [v,ix] = min(abs(ii-Ps));

        % Nearest prime number
        NP = Ps(ix);
        % Find nearest prime number based on node ID
        [d,e,n]=key_Generator(msk,NP);

        % Save keys
        memory(ii).dec = d;
        memory(ii).enc = e;
        memory(ii).n = n;

    end
    %%
    colors = {'r','g','b','m','y','c','k','r','g','b','m','y','c','r','g','b','m','y','c','k','r','g','b','m','y','c'};


    % Initalise energy
    E = Einit*ones(1,N);

    % Initilaise reading vector
    Avec = zeros(1,10000);
    Evec = zeros(1,10000);
    Dvec = zeros(1,10000);
    Vvec = zeros(1,10000);

    % Group of node which has not become CH for last 1/P round
    G = zeros(1,N);

    % First node died
    FND = 0;
    % Proposed


    % Start rounds
    for r11 = 1:10000
        speed1 = get(handles.checkbox1,'Value');

        % Calcualte number of alive nodes
        Avec(r11) = sum(E>0);

        % Dead nodes
        Dvec(r11) = N-Avec(r11);

        % if first dead node found then record it
        if Dvec(r11)>0 && FND==0
            FND = r11;
            fprintf('First node dead at round %d\n',r11);
            set(handles.edit10,'String',num2str(r11));
        end
        % Energy 
        Evec(r11) = sum(E);
        Vvec(r11) = var(E);


        fprintf('Round : %d, Alive nodes = %d, R Energy = %f\n',r11,Avec(r11),Evec(r11))
        set(handles.edit8,'String',num2str(r11))
        set(handles.edit7,'String',num2str(Avec(r11)));
        set(handles.edit9,'String',num2str(Evec(r11)));
        % If alive nodes are less then Nxp then stop simulation
        if Avec(r11)<N*p
            break;
        end


        if speed1 == 0
            %  --------------- Plot network ------------------------------
            % Plot the network
            plot(x,y,'k.');
            hold on

            % Plot basestation
            plot(BSx,BSy,'rs','markerfacecolor','g','markeredgecolor','k','markersize',12)
            text(BSx-3,BSy-5,'BS')

            axis([0 A+30 0 A+30]);



            % ------------------------------------------------------------
        end
    %% ======================== CH selection ==================================
        % Threshold caclulation
        T = (p/(1-p*mod(r11,round(1/p))))*E/Einit;

        % Get all the alive nodes
        Aix = find(E>0);

        % Get all the nodes who eligible to become CH
        Gi = find(G(Aix)<3);
        Ti = E(Gi);

        % For every eligible node generate random number and compare with
        % threshold
        Rand1 = rand(1,length(Gi));

        % Compare with threshold
        CHs = Aix(Gi(find(Rand1<Ti)));

        % If no cluseter head then 
        if isempty(CHs)

            % Sort G values
            [SortG sortix]  = sort(G(Aix));

            try
                % Take first N*p values
                CHs = Aix(sortix(1:round(N*p)));
            catch
                break;
            end
        elseif length(CHs)>10
            % Sort G values
            [SortG sortix]  = sort(G(Aix));

            % Take first N*p values
            CHs = Aix(sortix(1:round(N*p)));
        end

    %     
    % disp(length(CHs))
        % Set G counter for selected Cluster Head
        G(CHs) = 1/p;

        % GEt coordinates of cluster head
        CHx = x(CHs);
        CHy = y(CHs);

        % Divide nodes into clusters
        for ii = 1:N % For all nodes

            % GEt cordinates
            xt = x(ii);
            yt = y(ii);

            % Calculate distance to all cluster heads
            dists = sqrt((xt-CHx).^2+(yt-CHy).^2);

            % Find out nearest clusterhead
            [v ix] = min(dists);

            % update cluster ID
            Cid(ii) = ix;

            if speed1 == 0
                % Highlight with color
                plot(xt,yt,'o','markerfacecolor',colors{ix})
            end

        end

        if speed1 == 0
            % Plot cluster heads
            plot(CHx,CHy,'d','markersize',13);
        end
        % Peform communication of nodes to cluster head
        for ii = 1:N

            % GEt cordinates
            xt = x(ii);
            yt = y(ii);

            % Get cluster id
            ctemp = Cid(ii);

            % Calculate distance
            dist = sqrt((xt-CHx(ctemp))^2+(yt-CHy(ctemp))^2);

            % Calculate transmission energy
            Etx = calc_tx_energy_prop(dist,PL);

            E(ii) = E(ii)-Etx;

            % Calculate receiving energy
            Erx = calc_rx_energy_prop(PL);

            % REduce energy
            E(CHs(ctemp)) = E(CHs(ctemp)) - Erx;

            if E(ii)<=0
                continue;
            end

            % proposed
            % Sginature signing
            msg = ['NODE ID:' num2str(ii),', HELLO'];

            % Get encryption keys
            di = memory(ii).dec;
            ei = memory(ii).enc;
            ni = memory(ii).n ;

            % Generate Signature using homomorphic encryption
            [Sig1]= perform_encryption(msg,ei,ni);

            % Receiver side authenticate using decryption using decryption key
            rmsg = perform_decryption(Sig1,di,ni);

            if strcmpi(rmsg,msg)

                if speed1 == 0
                    % Communicate to respective cluster head
                    plot([xt CHx(ctemp)],[yt CHy(ctemp)],'c-');
                end 
            end
        end

        % Perform communication of cluster head with basestation
        for ii = 1:length(CHs)

            if speed1 == 0
                plot([CHx(ii) BSx],[CHy(ii) BSy],'r--','linewidth',2)
            end
            %  Calculate transmission energy
            dist = sqrt((CHx(ii)-BSx)^2+(CHy(ii)-BSy)^2);

            % Calculate number of cluster memeber
            cmi = find(Cid==ii);
            cm = sum(E(cmi)>0);

            % Calculate transmission energy
            Etx = calc_tx_energy_prop(dist,PL*cm);

            E(CHs(ii)) = E(CHs(ii))-Etx;


        end

        % Reduce the count for G
        G = G-1;

        % Make negative energy 0
        E(E<=0) = 0;
        if speed1 == 0
            % Highlight dead nodes
            Did = (E<=0);
            plot(x(Did),y(Did),'ks','markerfacecolor','k','markersize',5);
            hold off
            pause(0.001)
        end
        pause(0.0000001)
    end
    
    % Get the values
    Evec2 = Evec(1:r11);
    Avec2 = Avec(1:r11);
    Dvec2 = Dvec(1:r11);
    Vvec2 = Vvec(1:r11);
    E2vec2 = Evec2(1:end-1)-Evec2(2:end);
    E2vec2 = cumsum(E2vec2);
    FND2 = FND;

    save prop Evec2 Avec2 Dvec2 Vvec2 E2vec2 FND2
    % Plot the graphs
    figure;
    plot(Avec2)
    grid on
    xlabel('Rounds')
    ylabel('Alive nodes')
    title('Rounds vs. Alive nodes');

    figure;
    plot(Evec2)
    grid on
    xlabel('Rounds')
    ylabel('Residual Energy')
    title('Rounds vs. Residual energy');

    figure;
    plot(Dvec2)
    grid on
    xlabel('Rounds')
    ylabel('Dead nodes')
    title('Rounds vs. Dead nodes');

    figure;
    plot(Vvec2)
    grid on
    xlabel('Rounds')
    ylabel('Cost')
    title('Rounds vs. cost');

    figure;
    plot(E2vec2)
    grid on
    xlabel('Rounds')
    ylabel('Total energy consumption')
    title('Round vs. Total Energy consumption');
elseif get(handles.radiobutton2,'Value')==1
    
    msk = 971; % Any prime number

    % Generate memory for sensor nodes
    memory = struct;

   
    plot(x,y,'k.');
    hold on

    % Plot basestation
    plot(BSx,BSy,'rs','markerfacecolor','g','markeredgecolor','k','markersize',12)
    text(BSx-3,BSy-5,'BS')


    axis([0 A+30 0 A+30]);

    % Setup
    % Transmitt data to all nodes
    for ii = 1:N

        % -------------------------------------
        % Store in sensor node memory
        memory(ii).masterkey = msk;
        % Show communication line
        plot([x(ii) BSx],[y(ii) BSy],'r-');

    end
    title('Setup')
    pause(0.1)

    % Generate prime numbers
    Ps = primes(1000);
    Ps(1:2) = [];
    
    for ii = 1:N % For all nodes

        % Get master key
        msk = memory(ii).masterkey;

        % Find out nearest prime numbers
        [v,ix] = min(abs(ii-Ps));

        % Nearest prime number
        NP = Ps(ix);
        % Find nearest prime number based on node ID
        [d,e,n]=key_Generator(msk,NP);

        % Save keys
        memory(ii).dec = d;
        memory(ii).enc = e;
        memory(ii).n = n;

    end

    % Generate the signature offline
    for ii = 1:N
        % Proposed
            % Sginature signing
            msg = ['NODE ID:' num2str(ii),', HELLO'];

            % Get encryption keys
            di = memory(ii).dec;
            ei = memory(ii).enc;
            ni = memory(ii).n ;

            % Generate Signature using encryption
            [Sigoffline]= perform_encryption(msg,ei,ni);

            % Save offline signature
            memory(ii).Sigoffline=Sigoffline;
            memory(ii).msg = msg;
    end

    %%
    colors = {'r','g','b','m','y','c','k','r','g','b','m','y','c','r','g','b','m','y','c','k','r','g','b','m','y','c'};


    % Initalise energy
    E = Einit*ones(1,N);

    % Initilaise reading vector
    Avec = zeros(1,10000);
    Evec = zeros(1,10000);
    Dvec = zeros(1,10000);
    Vvec = zeros(1,10000);

    % Group of node which has not become CH for last 1/P round
    G = zeros(1,N);

    % First node died
    FND = 0;

    % Start rounds
    for r11 = 1:10000

        speed1 = get(handles.checkbox1,'Value');
        
        % Calcualte number of alive nodes
        Avec(r11) = sum(E>0);

        % Dead nodes
        Dvec(r11) = N-Avec(r11);

        % if first dead node found then record it
        if Dvec(r11)>0 && FND==0
            FND = r11;
            fprintf('First node dead at round %d\n',r11);
            set(handles.edit10,'String',num2str(r11));
        end
        % Energy 
        Evec(r11) = sum(E);
        Vvec(r11) = var(E);


        fprintf('Round : %d, Alive nodes = %d, R Energy = %f\n',r11,Avec(r11),Evec(r11))
        set(handles.edit8,'String',num2str(r11))
        set(handles.edit7,'String',num2str(Avec(r11)));
        set(handles.edit9,'String',num2str(Evec(r11)));
        % If alive nodes are less then Nxp then stop simulation
        if Avec(r11)<N*p
            break;
        end


        if speed1 == 0
            %  --------------- Plot network ------------------------------
            % Plot the network
            plot(x,y,'k.');
            hold on

            % Plot basestation
            plot(BSx,BSy,'rs','markerfacecolor','g','markeredgecolor','k','markersize',12)
            text(BSx-3,BSy-5,'BS')

            axis([0 A+30 0 A+30]);



        end
        % Threshold caclulation
        T = (p/(1-p*mod(r11,round(1/p))))*E/Einit;

        % Get all the alive nodes
        Aix = find(E>0);

        % Get all the nodes who eligible to become CH
        Gi = find(G(Aix)<3);
        Ti = E(Gi);

        Rand1 = rand(1,length(Gi));

        % Compare with threshold
        CHs = Aix(Gi(find(Rand1<Ti)));

        % If no cluseter head then 
        if isempty(CHs)

            % Sort G values
            [SortG sortix]  = sort(G(Aix));

            try
                % Take first N*p values
                CHs = Aix(sortix(1:round(N*p)));
            catch
                break;
            end
        elseif length(CHs)>10
            % Sort G values
            [SortG sortix]  = sort(G(Aix));

            % Take first N*p values
            CHs = Aix(sortix(1:round(N*p)));
        end

    %     
    % disp(length(CHs))
        % Set G counter for selected Cluster Head
        G(CHs) = 1/p;

        % GEt coordinates of cluster head
        CHx = x(CHs);
        CHy = y(CHs);

        % Divide nodes into clusters
        for ii = 1:N % For all nodes

            % GEt cordinates
            xt = x(ii);
            yt = y(ii);

            % Calculate distance to all cluster heads
            dists = sqrt((xt-CHx).^2+(yt-CHy).^2);

            % Find out nearest clusterhead
            [v ix] = min(dists);

            % update cluster ID
            Cid(ii) = ix;

            if speed1 == 0
                % Highlight with color
                plot(xt,yt,'o','markerfacecolor',colors{ix})
            end

        end

        if speed1 == 0
            % Plot cluster heads
            plot(CHx,CHy,'d','markersize',13);
        end
        % Peform communication of nodes to cluster head
        for ii = 1:N

            % GEt cordinates
            xt = x(ii);
            yt = y(ii);

            % Get cluster id
            ctemp = Cid(ii);

            % Calculate distance
            dist = sqrt((xt-CHx(ctemp))^2+(yt-CHy(ctemp))^2);

            % Calculate transmission energy
            Etx = calc_tx_energy_SLEACH(dist,PL);

            E(ii) = E(ii)-Etx;

            % Calculate receiving energy
            Erx = calc_rx_energy_SLEACH(PL);

            % REduce energy
            E(CHs(ctemp)) = E(CHs(ctemp)) - Erx;

            if E(ii)<=0
                continue;
            end
            % IBOOS

            % Load offline signature
            Sigoffline = memory(ii).Sigoffline;
            msg = memory(ii).msg;
            % Get decryption keys
            di = memory(ii).dec;
            ei = memory(ii).enc;
            ni = memory(ii).n ;


            % Receiver side authenticate using decryption using decryption key
            rmsg = perform_decryption(Sigoffline,di,ni);

            if strcmpi(rmsg,msg)
                if speed1 == 0
                    % Communicate to respective cluster head
                    plot([xt CHx(ctemp)],[yt CHy(ctemp)],'c-');
                end 
            end
        end

        % Perform communication of cluster head with basestation
        for ii = 1:length(CHs)

            if speed1 == 0
                plot([CHx(ii) BSx],[CHy(ii) BSy],'r--','linewidth',2)
            end
            %  Calculate transmission energy
            dist = sqrt((CHx(ii)-BSx)^2+(CHy(ii)-BSy)^2);

            % Calculate number of cluster memeber
            cmi = find(Cid==ii);
            cm = sum(E(cmi)>0);

            % Calculate transmission energy
            Etx = calc_tx_energy_SLEACH(dist,PL*cm);

            E(CHs(ii)) = E(CHs(ii))-Etx;


        end

        % Reduce the count for G
        G = G-1;

        % Make negative energy 0
        E(E<=0) = 0;
        if speed1 == 0
            % Highlight dead nodes
            Did = (E<=0);
            plot(x(Did),y(Did),'ks','markerfacecolor','k','markersize',5);
            hold off
            pause(0.001)
        end
    end

    % Get the values
    Evec3 = Evec(1:r11);
    Avec3 = Avec(1:r11);
    Dvec3 = Dvec(1:r11);
    Vvec3 = Vvec(1:r11);
    E2vec3 = Evec3(1:end-1)-Evec3(2:end);
    E2vec3 = cumsum(E2vec3);
    FND3 = FND;

    save SLEACH Evec3 Avec3 Dvec3 Vvec3 E2vec3 FND3
    % Plot the graphs
    figure;
    plot(Avec3)
    grid on
    xlabel('Rounds')
    ylabel('Alive nodes')
    title('Rounds vs. Alive nodes');

    figure;
    plot(Evec3)
    grid on
    xlabel('Rounds')
    ylabel('Residual Energy')
    title('Rounds vs. Residual energy');

    figure;
    plot(Dvec3)
    grid on
    xlabel('Rounds')
    ylabel('Dead nodes')
    title('Rounds vs. Dead nodes');

    figure;
    plot(Vvec3)
    grid on
    xlabel('Rounds')
    ylabel('Cost')
    title('Rounds vs. cost');

    figure;
    plot(E2vec3)
    grid on
    xlabel('Rounds')
    ylabel('Total energy consumption')
    title('Round vs. Total Energy consumption');
else
    colors = {'r','g','b','m','y','c','k','r','g','b','m','y','c','r','g','b','m','y','c','k','r','g','b','m','y','c'};


% Initalise energy
E = Einit*ones(1,N);

% Initilaise reading vector
Avec = zeros(1,10000);
Evec = zeros(1,10000);
Dvec = zeros(1,10000);
Vvec = zeros(1,10000);


% Group of node which has not become CH for last 1/P round
G = zeros(1,N);

% First node died
FND = 0;

% Start rounds
for r11 = 1:10000
    
    speed1 = get(handles.checkbox1,'Value');
    
    % Calcualte number of alive nodes
    Avec(r11) = sum(E>0);
    
    % Dead nodes
    Dvec(r11) = N-Avec(r11);
    
    % if first dead node found then record it
    if Dvec(r11)>0 && FND==0
        FND = r11;
        fprintf('First node dead at round %d\n',r11);
        set(handles.edit10,'String',num2str(r11));
    end
    % Energy 
    Evec(r11) = sum(E);
    Vvec(r11) = var(E);
    fprintf('Round : %d, Alive nodes = %d, R Energy = %f\n',r11,Avec(r11),Evec(r11))
    set(handles.edit8,'String',num2str(r11))
    set(handles.edit7,'String',num2str(Avec(r11)));
    set(handles.edit9,'String',num2str(Evec(r11)));
    % If alive nodes are less then Nxp then stop simulation
    if Avec(r11)<N*p
        break;
    end
    if speed1 == 0
        plot(x,y,'k.');
        hold on

        % Plot basestation
        plot(BSx,BSy,'rs','markerfacecolor','g','markeredgecolor','k','markersize',12)
        text(BSx-3,BSy-5,'BS')
       
        axis([0 A+30 0 A+30]);
    end

    % Threshold caclulation
    T = (p/(1-p*mod(r11,round(1/p))))*E/Einit;

    % Get all the alive nodes
    Aix = find(E>0);
    
    % Get all the nodes who eligible to become CH
    Gi = find(G(Aix)<3);
    Ti = E(Gi);
    Rand1 = rand(1,length(Gi));
    
    % Compare with threshold
    CHs = Aix(Gi(find(Rand1<Ti)));
    
    % If no cluseter head then 
    if isempty(CHs)
        
        % Sort G values
        [SortG sortix]  = sort(G(Aix));
        
        try
            % Take first N*p values
            CHs = Aix(sortix(1:round(N*p)));
        catch
            break;
        end
    elseif length(CHs)>10
        % Sort G values
        [SortG sortix]  = sort(G(Aix));
        
        % Take first N*p values
        CHs = Aix(sortix(1:round(N*p)));
    end
    

    G(CHs) = 1/p;
    
    % GEt coordinates of cluster head
    CHx = x(CHs);
    CHy = y(CHs);
    
    % Divide nodes into clusters
    for ii = 1:N % For all nodes
        
        % GEt cordinates
        xt = x(ii);
        yt = y(ii);
        
        % Calculate distance to all cluster heads
        dists = sqrt((xt-CHx).^2+(yt-CHy).^2);
        
        % Find out nearest clusterhead
        [v ix] = min(dists);
        
        % update cluster ID
        Cid(ii) = ix;
        
        if speed1 == 0
            % Highlight with color
            plot(xt,yt,'o','markerfacecolor',colors{ix})
        end
 
    end
    
    if speed1 == 0
        % Plot cluster heads
        plot(CHx,CHy,'d','markersize',13);
    end
    % Peform communication of nodes to cluster head
    for ii = 1:N
        
        % GEt cordinates
        xt = x(ii);
        yt = y(ii);
        
        % Get cluster id
        ctemp = Cid(ii);
        
        % Calculate distance
        dist = sqrt((xt-CHx(ctemp))^2+(yt-CHy(ctemp))^2);
        
        % Calculate transmission energy
        Etx = calc_tx_energy(dist,PL);
        
        E(ii) = E(ii)-Etx;
        
        % Calculate receiving energy
        Erx = calc_rx_energy(PL);
        
        % REduce energy
        E(CHs(ctemp)) = E(CHs(ctemp)) - Erx;
        
        if E(ii)<=0
            continue;
        end
        if speed1 == 0
            % Communicate to respective cluster head
            plot([xt CHx(ctemp)],[yt CHy(ctemp)],'c-');
        end 
    end

    % Perform communication of cluster head with basestation
    for ii = 1:length(CHs)
        
        if speed1 == 0
            plot([CHx(ii) BSx],[CHy(ii) BSy],'r--','linewidth',2)
        end
        %  Calculate transmission energy
        dist = sqrt((CHx(ii)-BSx)^2+(CHy(ii)-BSy)^2);
        
        % Calculate number of cluster memeber
        cmi = find(Cid==ii);
        cm = sum(E(cmi)>0);
        
        % Calculate transmission energy
        Etx = calc_tx_energy(dist,PL*cm);
        
        E(CHs(ii)) = E(CHs(ii))-Etx;     
    end

    % Reduce the count for G
    G = G-1;
    
    % Make negative energy 0
    E(E<=0) = 0;
    if speed1 == 0
        % Highlight dead nodes
        Did = (E<=0);
        plot(x(Did),y(Did),'ks','markerfacecolor','k','markersize',5);
        hold off
        pause(0.001)
    end
end

% Get the values
Evec1 = Evec(1:r11);
Avec1 = Avec(1:r11);
Dvec1 = Dvec(1:r11);
Vvec1 = Vvec(1:r11);
E2vec1 = Evec1(1:end-1)-Evec1(2:end);
E2vec1 = cumsum(E2vec1);
FND1 = FND;

save SecLEACH Evec1 Avec1 Dvec1 Vvec1 E2vec1 FND1
% Plot the graphs
figure;
plot(Avec1)
grid on
xlabel('Rounds')
ylabel('Alive nodes')
title('Rounds vs. Alive nodes');

figure;
plot(Evec1)
grid on
xlabel('Rounds')
ylabel('Residual Energy')
title('Rounds vs. Residual energy');

figure;
plot(Dvec1)
grid on
xlabel('Rounds')
ylabel('Dead nodes')
title('Rounds vs. Dead nodes');

figure;
plot(Vvec1)
grid on
xlabel('Rounds')
ylabel('Cost')
title('Rounds vs. cost');

figure;
plot(E2vec1)
grid on
xlabel('Rounds')
ylabel('Total energy consumption')
title('Round vs. Total Energy consumption');
end

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)

function edit7_Callback(hObject, eventdata, handles)

function edit7_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit8_Callback(hObject, eventdata, handles)

function edit8_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit9_Callback(hObject, eventdata, handles)
function edit9_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
clc;
clear all;
close all;
%%
load SecLEACH
load prop
load SLEACH
figure;
bar(1,FND1)
hold on
bar(2,FND2,'r')
bar(3,FND3,'g')
legend('SecLEACH','Proposed','S-LEACH')
xlabel('Protocol');
ylabel('Memory Dependency');

% Plot the graphs
figure;
plot(Avec1,'r-')
hold on
plot(Avec2,'b-')
plot(Avec3,'g-')
grid on
xlabel('Rounds')
ylabel('Alive nodes')
title('Rounds vs. Alive nodes');
legend('Prop','SecLEACH','S-LEACH')

figure;
plot(Evec1,'r-')
hold on
plot(Evec2,'b-')
plot(Evec3,'g-')
grid on
xlabel('Rounds')
ylabel('Residual Energy')
title('Rounds vs. Residual energy');
legend('Prop','SecLEACH','S-LEACH')

figure;
plot(Vvec1)
hold on
plot(Vvec2,'r-')
plot(Vvec3,'g-')
grid on
xlabel('Rounds')
ylabel('Cost')
title('Rounds vs. cost');
legend('SecLeach','Prop','S-LEACH')


figure;
plot(E2vec1)
hold on
plot(E2vec2,'r-')
plot(E2vec3,'g-')
grid on
xlabel('Rounds')
ylabel('Total energy consumption')
title('Round vs. Total Energy consumption');
legend('SecLeach','Prop','S-LEACH')

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
close(handles.figure1)
function edit10_Callback(hObject, eventdata, handles)

function edit10_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
