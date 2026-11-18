#!/usr/bin/env python3
"""hop_depth.py file.gr... -- hop eccentricity (BFS) of node 1 and hop depth (max, mean) of the
shortest-path tree from node 1 (Dijkstra, ties broken by first arrival). Used for the analysis of
first-order methods on shortest paths (results/depth/)."""
import heapq, sys, collections
def load(p):
    adj=None
    for l in open(p):
        if l[0]=='p': n=int(l.split()[2]); adj=[[] for _ in range(n+1)]
        elif l[0]=='a':
            _,u,v,w=l.split(); adj[int(u)].append((int(v),int(w)))
    return adj
for p in sys.argv[1:]:
    adj=load(p); n=len(adj)-1; m=sum(map(len,adj))
    # BFS eccentricity of node 1
    d=[-1]*(n+1); d[1]=0; q=collections.deque([1])
    while q:
        u=q.popleft()
        for v,_ in adj[u]:
            if d[v]<0: d[v]=d[u]+1; q.append(v)
    ecc=max(d)
    # Dijkstra; hop depth along the shortest-path tree
    D=[float('inf')]*(n+1); H=[0]*(n+1); D[1]=0; pq=[(0,1)]
    while pq:
        du,u=heapq.heappop(pq)
        if du>D[u]: continue
        for v,w in adj[u]:
            if du+w<D[v]: D[v]=du+w; H[v]=H[u]+1; heapq.heappush(pq,(D[v],v))
    reach=[h for h,x in zip(H[1:],D[1:]) if x<float('inf')]
    print(f"{p.split('/')[-1]:28s} n={n:7d} m={m:8d} bfs_ecc={ecc:5d} spt_depth={max(reach):5d} spt_mean={sum(reach)/len(reach):7.1f}")
