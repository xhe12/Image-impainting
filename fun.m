function [FDx,FDy,BDx,BDy,CDx,CDy,avgx,avgy,Ac,As] = fun

FDx = @(U) ForwardDx(U);
FDy = @(U) ForwardDy(U);
BDx = @(U)BackwardDx(U);
BDy = @(U)BackwardDy(U);
CDx = @(U)CentralDx(U);
CDy = @(U)CentralDy(U);
avgx = @(U)ax(U);
avgy = @(U)ay(U);
Ac = @(U)Avgc(U);
As = @(U)Avgs(U);

    function Dux = ForwardDx(U)
    Dux = [U(:,2:end)-U(:,1:end-1), U(:,1) - U(:,end)];
    end

    function Duy = ForwardDy(U)
    Duy = [U(2:end,:)-U(1:end-1,:); U(1,:) - U(end,:)];
    end

    function BDux = BackwardDx(U)
        BDux = [U(:,1) - U(:,end), U(:,2:end)-U(:,1:end-1)];
    end

    function BDuy = BackwardDy(U)
        BDuy = [U(1,:) - U(end,:); U(2:end,:)-U(1:end-1,:)];
    end

    function CDux = CentralDx(U)
        CDux = [U(:,2)-U(:,end), U(:,3:end)-U(:,1:end-2), U(:,1)-U(:,end-1)];
        CDux = CDux/2;
    end

    function CDuy = CentralDy(U)
        CDuy = [U(2,:)-U(end,:); U(3:end,:)-U(1:end-2,:); U(1,:)-U(end-1,:)];
        CDuy = CDuy/2;
    end


    function avg = ax(U)
        %Average in the x direction u(i,j) = (u(i,j)+u(i,j+1))/2;
        avg = [(U(:,1:end-1)+U(:,2:end))/2 (U(:,end)+U(:,1))/2];
    end

    function avg = ay(U)
        %Average in the y direction u(i,j) = (u(i,j)+u(i+1,j))/2;
        avg = [(U(1:end-1,:)+U(2:end,:))/2; (U(end,:)+U(1,:))/2];
    end

    function avg = Avgs(U)
        [Y,X] = size(U);
        %Enlarge U;
        U = repmat(U,3);
        U = U(Y:2*Y+1,X:2*X+1);
        % Transfer the value on the circle nodes to the square nodes;
        avg = (U(3:Y+2,2:X+1)+U(3:Y+2,1:X)+U(2:Y+1,2:X+1)+U(2:Y+1,1:X))/4;
    end

    function avg = Avgc(U)
        [Y,X] = size(U);
        %Enlarge U;
        U = repmat(U,3);
        U = U(Y:2*Y+1,X:2*X+1);
        % Transfer the value on the square nodes to the circle nodes;
        avg = (U(2:Y+1,3:X+2)+U(2:Y+1,2:X+1)+U(1:Y,3:X+2)+U(1:Y,2:X+1))/4;
    end
        
end

