function [m_newM,m_newS,StD,V]=RunningStat(iter,x,m_oldM,m_oldS)
% DEBUG
%iter = irep-nburn;
%x = deltaB.^2;
%RunningStat(irep-nburn,deltaB.^2,[],[])


if iter==1
    m_newM = x;
    m_newS = 0.0;            
else
    m_newM = m_oldM + (x - m_oldM)./iter;
    m_newS = m_oldS + (x - m_oldM).*(x - m_newM);
end
V = m_newS./(iter-1);
StD = sqrt(V);