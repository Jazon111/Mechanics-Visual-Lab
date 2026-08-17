
const NODES = { A:{id:'A',x:0,y:0}, B:{id:'B',x:6,y:0}, C:{id:'C',x:0,y:4}, D:{id:'D',x:6,y:4}, E:{id:'E',x:8,y:4} };
const MEMBERS = [
  { id:'AC', nodeA:'A', nodeB:'C', L:4, angle:Math.PI/2 },
  { id:'CD', nodeA:'C', nodeB:'D', L:6, angle:0 },
  { id:'DE', nodeA:'D', nodeB:'E', L:2, angle:0 },
  { id:'BD', nodeA:'B', nodeB:'D', L:4, angle:Math.PI/2 }
];
const SUPPORTS = { A:[true,true,false], B:[false,true,false], C:[false,false,false], D:[false,false,false], E:[false,false,false] };
const NODE_IDS = ['A','B','C','D','E'];
const NODES_ALL = NODE_IDS.length;
const DOF_TOTAL = NODES_ALL * 3;

const dofMap = {};
NODE_IDS.forEach((id, i) => { dofMap[id] = [i*3, i*3+1, i*3+2]; });

// 构建约束索引
const constrainedDOF = [];
NODE_IDS.forEach(id => {
  SUPPORTS[id].forEach((c, i) => { if(c) constrainedDOF.push(dofMap[id][i]); });
});
const freeDOF = [];
for(let i=0;i<DOF_TOTAL;i++){ if(!constrainedDOF.includes(i)) freeDOF.push(i); }

// 材料参数
const E = 2e8;
const A_sec = 0.01;
const I_sec = 0.0001;

// 初始荷载
let loads = [
  { id:1, type:'udl', member:'CD', value:10, pos1:0, pos2:6 },
  { id:2, type:'point', member:'DE', pos:2, value:9, direction:'y' },
  { id:3, type:'point', member:'BD', pos:2, value:-6, direction:'x' }
];

// Canvas
const canvas = document.getElementById('scene');
const ctx = canvas.getContext('2d');
let W=0, H=0, SCALE=55, ox=0, oy=0;
let showM=true, showV=false, showN=false, showFrame=true;
let solved=false, solution=null;

function resize(){
  const rect = canvas.getBoundingClientRect();
  const dpr = Math.min(window.devicePixelRatio||1,2);
  canvas.width = Math.floor(rect.width)*dpr;
  canvas.height = Math.floor(rect.height)*dpr;
  ctx.setTransform(dpr,0,0,dpr,0,0);
  W=rect.width; H=rect.height;
  // 适应实际刚架尺寸: x范围0-8, y范围0-4
  const padL=50, padR=30, padT=40, padB=60;
  SCALE = Math.min((W-padL-padR)/8, (H-padT-padB)/4);
  SCALE = Math.max(8, Math.min(80, SCALE));
  ox = padL; oy = H-padB;
  render();
}
window.addEventListener('resize', resize);

function toS(x,y){ return {x:ox+x*SCALE, y:oy-y*SCALE}; }

// ========== 矩阵位移法 ==========
function solve(){
  // 1. 组装整体刚度矩阵 (DOF_TOTAL x DOF_TOTAL)
  const K = Array.from({length:DOF_TOTAL}, () => Array(DOF_TOTAL).fill(0));
  
  MEMBERS.forEach(mem => {
    const Ke = elementK(mem);
    const dofs = [...dofMap[mem.nodeA], ...dofMap[mem.nodeB]];
    for(let i=0;i<6;i++) for(let j=0;j<6;j++){
      K[dofs[i]][dofs[j]] += Ke[i][j];
    }
  });
  
  // 2. 计算等效节点荷载
  const F = computeEquivalentLoads();
  
  // 3. 提取自由自由度的子矩阵
  const nFree = freeDOF.length;
  const Kff = Array.from({length:nFree}, () => Array(nFree).fill(0));
  const Ff = Array(nFree).fill(0);
  for(let i=0;i<nFree;i++){
    Ff[i] = F[freeDOF[i]];
    for(let j=0;j<nFree;j++){
      Kff[i][j] = K[freeDOF[i]][freeDOF[j]];
    }
  }
  
  // 4. 求解 Kff * Uf = Ff
  const Uf = gaussSolve(Kff, Ff);
  
  // 5. 组装完整位移向量
  const U = Array(DOF_TOTAL).fill(0);
  for(let i=0;i<nFree;i++) U[freeDOF[i]] = Uf[i];
  
  // 6. 计算反力 (约束自由度)
  const reactions = {};
  NODE_IDS.forEach(id => {
    const sup = SUPPORTS[id];
    const dofs = dofMap[id];
    reactions[id] = { fx:0, fy:0, m:0 };
    if(sup[0]){
      let r=0; for(let j=0;j<DOF_TOTAL;j++) r+=K[dofs[0]][j]*U[j];
      reactions[id].fx = r - F[dofs[0]];
    }
    if(sup[1]){
      let r=0; for(let j=0;j<DOF_TOTAL;j++) r+=K[dofs[1]][j]*U[j];
      reactions[id].fy = r - F[dofs[1]];
    }
    if(sup[2]){
      let r=0; for(let j=0;j<DOF_TOTAL;j++) r+=K[dofs[2]][j]*U[j];
      reactions[id].m = r - F[dofs[2]];
    }
  });
  
  // 7. 计算各杆件内力
  const memberForces = computeInternalForces(U);
  
  return { U, reactions, memberForces };
}

function elementK(mem){
  const L=mem.L, angle=mem.angle;
  const c=Math.cos(angle), s=Math.sin(angle);
  const eal=E*A_sec/L, eil=E*I_sec/L;
  const eil2=E*I_sec/(L*L), eil3=E*I_sec/(L*L*L);
  
  const kl = [
    [eal,0,0,-eal,0,0],
    [0,12*eil3,6*eil2,0,-12*eil3,6*eil2],
    [0,6*eil2,4*eil,0,-6*eil2,2*eil],
    [-eal,0,0,eal,0,0],
    [0,-12*eil3,-6*eil2,0,12*eil3,-6*eil2],
    [0,6*eil2,2*eil,0,-6*eil2,4*eil]
  ];
  
  const T = [[c,s,0,0,0,0],[-s,c,0,0,0,0],[0,0,1,0,0,0],[0,0,0,c,s,0],[0,0,0,-s,c,0],[0,0,0,0,0,1]];
  const Tt = transpose(T);
  return matMul(matMul(Tt,kl),T);
}

function matMul(A,B){
  const n=A.length,m=B[0].length,p=B.length,C=Array.from({length:n},()=>Array(m).fill(0));
  for(let i=0;i<n;i++) for(let j=0;j<m;j++) for(let k=0;k<p;k++) C[i][j]+=A[i][k]*B[k][j];
  return C;
}
function transpose(A){ return A[0].map((_,i)=>A.map(r=>r[i])); }

function computeEquivalentLoads(){
  const F = Array(DOF_TOTAL).fill(0);
  
  loads.forEach(ld => {
    const mem = MEMBERS.find(m=>m.id===ld.member);
    if(!mem) return;
    const L=mem.L, angle=mem.angle;
    const c=Math.cos(angle), s=Math.sin(angle);
    const dA=dofMap[mem.nodeA], dB=dofMap[mem.nodeB];
    
    if(ld.type==='udl'){
      const w=ld.value, a=ld.pos1, b=ld.pos2, Lspan=b-a;
      const FyA=-w*Lspan/2, FyB=-w*Lspan/2;
      const mA=-w*Lspan*Lspan/12, mB=w*Lspan*Lspan/12;
      // 局部坐标力 (反向): N=0, Vy=-FyA, M=-mA
      const fA=[0,-FyA,-mA], fB=[0,-FyB,-mB];
      // 旋转到整体坐标
      const fAG = [c*fA[0]-s*fA[1], s*fA[0]+c*fA[1], fA[2]];
      const fBG = [c*fB[0]-s*fB[1], s*fB[0]+c*fB[1], fB[2]];
      [dA,dB].forEach((dofs,idx)=>{
        const f = idx===0?fAG:fBG;
        for(let i=0;i<3;i++) F[dofs[i]] -= f[i];
      });
    } else if(ld.type==='point'){
      const P=ld.value, a=ld.pos, b=L-a;
      // 全局坐标方向转局部坐标
      // direction='x': 水平力(+向右) → 局部: N=c*P, V=-s*P
      // direction='y': 竖力(+向下) → 全局(0,-P) → 局部: N=-s*P, V=-c*P
      let Nloc, Vloc;
      if(ld.direction==='y'){ Nloc=-s*P; Vloc=-c*P; }
      else { Nloc=c*P; Vloc=-s*P; }
      // 轴力部分 (无固端弯矩)
      const NA=Nloc*b/L, NB=Nloc*a/L;
      // 剪力部分
      const VA=Vloc*b/L, VB=Vloc*a/L;
      const mA_fem=-Vloc*a*b*b/(L*L), mB_fem=Vloc*a*a*b/(L*L);
      // 局部杆端力 (反向)
      const fA=[-NA,-VA,-mA_fem], fB=[-NB,-VB,-mB_fem];
      const fAG=[c*fA[0]-s*fA[1],s*fA[0]+c*fA[1],fA[2]];
      const fBG=[c*fB[0]-s*fB[1],s*fB[0]+c*fB[1],fB[2]];
      for(let i=0;i<3;i++){ F[dA[i]]-=fAG[i]; F[dB[i]]-=fBG[i]; }
    } else if(ld.type==='moment'){
      const m=ld.value, a=ld.pos, b=L-a;
      const mA=m*b*(3*a+b)/(L*L), mB=-m*a*(a+3*b)/(L*L);
      const VyA=-6*m*a*b/(L*L*L), VyB=6*m*a*b/(L*L*L);
      const fA=[0,-VyA,-mA], fB=[0,-VyB,-mB];
      const fAG=[c*fA[0]-s*fA[1],s*fA[0]+c*fA[1],fA[2]];
      const fBG=[c*fB[0]-s*fB[1],s*fB[0]+c*fB[1],fB[2]];
      for(let i=0;i<3;i++){ F[dA[i]]-=fAG[i]; F[dB[i]]-=fBG[i]; }
    }
  });
  
  return F;
}

function computeInternalForces(U){
  const results = {};
  MEMBERS.forEach(mem => {
    const L=mem.L, angle=mem.angle;
    const c=Math.cos(angle), s=Math.sin(angle);
    const dofsA=dofMap[mem.nodeA], dofsB=dofMap[mem.nodeB];
    
    // 节点位移 (整体坐标)
    const uA=[U[dofsA[0]],U[dofsA[1]],U[dofsA[2]]];
    const uB=[U[dofsB[0]],U[dofsB[1]],U[dofsB[2]]];
    
    // 转到局部坐标
    const uAL=c*uA[0]+s*uA[1], vAL=-s*uA[0]+c*uA[1], θAL=uA[2];
    const uBL=c*uB[0]+s*uB[1], vBL=-s*uB[0]+c*uB[1], θBL=uB[2];
    
    // 局部杆端力 (由位移直接计算, 已包含荷载效应)
    const eal=E*A_sec/L, eil=E*I_sec/L, eil2=E*I_sec/(L*L), eil3=E*I_sec/(L*L*L);
    const NA = eal*(uBL-uAL);
    const VyA = 12*eil3*vAL+6*eil2*θAL-12*eil3*vBL+6*eil2*θBL;
    const MA = 6*eil2*vAL+4*eil*θAL-6*eil2*vBL+2*eil*θBL;
    const NB = -NA;
    const VyB = -12*eil3*vAL-6*eil2*θAL+12*eil3*vBL-6*eil2*θBL;
    const MB = 6*eil2*vAL+2*eil*θAL-6*eil2*vBL+4*eil*θBL;
    
    // 计算沿杆的内力分布
    const pts=80, xs=[], Ms=[], Vs=[], Ns=[];
    for(let i=0;i<=pts;i++){
      const x=(i/pts)*L; xs.push(x);
      let m=MA+VyA*x, v=VyA, n=NA;
      // 叠加荷载效应
      loads.filter(l=>l.member===mem.id).forEach(ld=>{
        if(ld.type==='udl'){
          const a=ld.pos1, Lspan=ld.pos2-ld.pos1, w=ld.value;
          if(x>=a){
            const dx=Math.min(x-a, Lspan);
            m -= w*dx*dx/2; v -= w*dx;
          }
        } else if(ld.type==='point'){
          if(x>ld.pos){
            let Nloc, Vloc;
            if(ld.direction==='y'){ Nloc=-s*ld.value; Vloc=-c*ld.value; }
            else { Nloc=c*ld.value; Vloc=-s*ld.value; }
            m -= Vloc*(x-ld.pos);
            v -= Vloc;
            n -= Nloc;
          }
        }
      });
      Ms.push(m); Vs.push(v); Ns.push(n);
    }
    
    results[mem.id] = { M:Ms, V:Vs, N:Ns, x:xs, L, angle, 
      start:{N:NA,V:VyA,M:MA}, end:{N:NB,V:VyB,M:MB} };
  });
  return results;
}

function gaussSolve(A,b){
  const n=A.length, M=A.map((r,i)=>[...r,b[i]]);
  for(let i=0;i<n;i++){
    let mr=i;
    for(let k=i+1;k<n;k++) if(Math.abs(M[k][i])>Math.abs(M[mr][i])) mr=k;
    [M[i],M[mr]]=[M[mr],M[i]];
    for(let k=i+1;k<n;k++){
      const f=M[k][i]/M[i][i];
      for(let j=i;j<=n;j++) M[k][j]-=f*M[i][j];
    }
  }
  const x=Array(n).fill(0);
  for(let i=n-1;i>=0;i--){
    x[i]=M[i][n];
    for(let j=i+1;j<n;j++) x[i]-=M[i][j]*x[j];
    x[i]/=M[i][i];
  }
  return x;
}

// ========== Canvas 渲染 ==========
function render(){
  ctx.clearRect(0,0,W,H);
  drawGrid();
  if(showFrame) drawStructure();
  if(!solved) drawLoads();
  if(solved && solution){
    if(showM) drawDiagram(solution.memberForces,'M','#3fb950');
    if(showV) drawDiagram(solution.memberForces,'V','#ffd700');
    if(showN) drawDiagram(solution.memberForces,'N','#a371f7');
    drawReactions(solution.reactions);
  }
}

function drawGrid(){
  ctx.save();
  ctx.strokeStyle='rgba(88,166,255,0.05)';
  ctx.lineWidth=1;
  for(let x=ox%SCALE;x<W;x+=SCALE){ ctx.beginPath(); ctx.moveTo(x,0); ctx.lineTo(x,H); ctx.stroke(); }
  for(let y=oy%SCALE;y<H;y+=SCALE){ ctx.beginPath(); ctx.moveTo(0,y); ctx.lineTo(W,y); ctx.stroke(); }
  ctx.restore();
}

function drawStructure(){
  // 杆件
  MEMBERS.forEach(mem => {
    const sa=toS(NODES[mem.nodeA].x, NODES[mem.nodeA].y);
    const sb=toS(NODES[mem.nodeB].x, NODES[mem.nodeB].y);
    ctx.save();
    ctx.strokeStyle='#30363d'; ctx.lineWidth=14; ctx.lineCap='round';
    ctx.beginPath(); ctx.moveTo(sa.x,sa.y); ctx.lineTo(sb.x,sb.y); ctx.stroke();
    ctx.strokeStyle='#484f58'; ctx.lineWidth=10;
    ctx.beginPath(); ctx.moveTo(sa.x,sa.y); ctx.lineTo(sb.x,sb.y); ctx.stroke();
    ctx.strokeStyle='#6e7681'; ctx.lineWidth=6;
    ctx.shadowColor='rgba(88,166,255,0.3)'; ctx.shadowBlur=3;
    ctx.beginPath(); ctx.moveTo(sa.x,sa.y); ctx.lineTo(sb.x,sb.y); ctx.stroke();
    ctx.restore();
  });
  
  // 节点
  NODE_IDS.forEach(id => {
    const s=toS(NODES[id].x, NODES[id].y);
    ctx.save();
    ctx.fillStyle='#0d1117'; ctx.strokeStyle='#58a6ff'; ctx.lineWidth=2;
    ctx.beginPath(); ctx.arc(s.x,s.y,7,0,Math.PI*2); ctx.fill(); ctx.stroke();
    ctx.fillStyle='#58a6ff';
    ctx.beginPath(); ctx.arc(s.x,s.y,3,0,Math.PI*2); ctx.fill();
    ctx.font='bold 13px "Segoe UI"'; ctx.fillStyle='#c9d1d9'; ctx.textAlign='center';
    ctx.fillText(id, s.x+(s.x<W/2?-18:18), s.y+(s.y<H/2?-12:14));
    ctx.restore();
  });
  
  // 支座
  const aS=toS(0,0), bS=toS(6,0);
  ctx.save(); ctx.fillStyle='#58a6ff';
  // A 铰
  ctx.beginPath(); ctx.moveTo(aS.x,aS.y); ctx.lineTo(aS.x-15,aS.y+22); ctx.lineTo(aS.x+15,aS.y+22); ctx.closePath(); ctx.fill();
  ctx.strokeStyle='#58a6ff'; ctx.lineWidth=2;
  ctx.beginPath(); ctx.moveTo(aS.x-22,aS.y+22); ctx.lineTo(aS.x+22,aS.y+22); ctx.stroke();
  for(let i=-3;i<=3;i++){ ctx.beginPath(); ctx.moveTo(aS.x+i*6-3,aS.y+22); ctx.lineTo(aS.x+i*6+3,aS.y+28); ctx.stroke(); }
  // B 可动铰
  ctx.fillStyle='#58a6ff';
  ctx.beginPath(); ctx.moveTo(bS.x,bS.y); ctx.lineTo(bS.x-15,bS.y+18); ctx.lineTo(bS.x+15,bS.y+18); ctx.closePath(); ctx.fill();
  ctx.beginPath(); ctx.arc(bS.x-10,bS.y+24,4,0,Math.PI*2); ctx.arc(bS.x+10,bS.y+24,4,0,Math.PI*2); ctx.fill();
  ctx.strokeStyle='#58a6ff';
  ctx.beginPath(); ctx.moveTo(bS.x-22,bS.y+30); ctx.lineTo(bS.x+22,bS.y+30); ctx.stroke();
  for(let i=-3;i<=3;i++){ ctx.beginPath(); ctx.moveTo(bS.x+i*6-3,bS.y+30); ctx.lineTo(bS.x+i*6+3,bS.y+36); ctx.stroke(); }
  ctx.restore();
  
  // 尺寸标注
  ctx.save(); ctx.font='bold 11px "Consolas",monospace'; ctx.fillStyle='#8b949e'; ctx.textAlign='center';
  MEMBERS.forEach(mem => {
    const sa=toS(NODES[mem.nodeA].x, NODES[mem.nodeA].y);
    const sb=toS(NODES[mem.nodeB].x, NODES[mem.nodeB].y);
    const mx=(sa.x+sb.x)/2, my=(sa.y+sb.y)/2;
    if(mem.angle===0){ ctx.fillText(mem.L+'m', mx, my+18); }
    else { ctx.fillText(mem.L+'m', mx+16, my+4); }
  });
  // 跨度
  const aT=toS(0,0), bT=toS(6,0);
  ctx.fillText('6m', (aT.x+bT.x)/2, aT.y+48);
  const dT=toS(6,4), eT=toS(8,4);
  ctx.fillText('2m', (dT.x+eT.x)/2, dT.y-14);
  ctx.restore();
}

function drawLoads(){
  loads.forEach(ld => {
    const mem=MEMBERS.find(m=>m.id===ld.member); if(!mem) return;
    const sa=toS(NODES[mem.nodeA].x, NODES[mem.nodeA].y);
    const sb=toS(NODES[mem.nodeB].x, NODES[mem.nodeB].y);
    const angle=Math.atan2(sb.y-sa.y, sb.x-sa.x);
    
    if(ld.type==='udl'){
      const p1=ld.pos1/mem.L, p2=ld.pos2/mem.L;
      const s1x=sa.x+(sb.x-sa.x)*p1, s1y=sa.y+(sb.y-sa.y)*p1;
      const s2x=sa.x+(sb.x-sa.x)*p2, s2y=sa.y+(sb.y-sa.y)*p2;
      const len=Math.sqrt((s2x-s1x)**2+(s2y-s1y)**2);
      const cnt=Math.max(3,Math.floor(len/18));
      const sp=len/cnt;
      // 垂直方向 (朝下为正)
      const perpAngle = angle + Math.PI/2;
      const off = 30;
      ctx.save();
      ctx.strokeStyle='#ffd700'; ctx.lineWidth=2;
      ctx.shadowColor='#ffd700'; ctx.shadowBlur=4;
      // 顶部线
      const offS = {x:Math.cos(perpAngle+Math.PI)*off, y:Math.sin(perpAngle+Math.PI)*off};
      const offE = {x:Math.cos(perpAngle+Math.PI)*off, y:Math.sin(perpAngle+Math.PI)*off};
      ctx.beginPath(); ctx.moveTo(s1x+offS.x,s1y+offS.y); ctx.lineTo(s2x+offE.x,s2y+offE.y); ctx.stroke();
      // 箭头
      for(let i=0;i<=cnt;i++){
        const t=p1+(p2-p1)*i/cnt;
        const px=sa.x+(sb.x-sa.x)*t, py=sa.y+(sb.y-sa.y)*t;
        const tipX=px+Math.cos(perpAngle+Math.PI)*off, tipY=py+Math.sin(perpAngle+Math.PI)*off;
        const baseX=px+Math.cos(perpAngle)*12, baseY=py+Math.sin(perpAngle)*12;
        ctx.beginPath(); ctx.moveTo(tipX,tipY); ctx.lineTo(baseX,baseY); ctx.stroke();
      }
      ctx.shadowBlur=0;
      ctx.font='bold 11px "Consolas",monospace'; ctx.fillStyle='#ffd700'; ctx.textAlign='center';
      const midX=(s1x+s2x)/2+offS.x*0.3, midY=(s1y+s2y)/2+offS.y*0.3;
      ctx.fillText(ld.value+' kN/m', midX, midY);
      ctx.restore();
    } else if(ld.type==='point'){
      const t=ld.pos/mem.L;
      const px=sa.x+(sb.x-sa.x)*t, py=sa.y+(sb.y-sa.y)*t;
      let dirAng;
      if(ld.direction==='y'){ dirAng=Math.PI/2; }
      else { dirAng=0; }
      const sign = ld.value>=0?1:-1;
      const tipX=px+Math.cos(dirAng)*24*sign, tipY=py+Math.sin(dirAng)*24*sign;
      ctx.save();
      ctx.strokeStyle='#f85149'; ctx.fillStyle='#f85149'; ctx.lineWidth=2.5;
      ctx.shadowColor='#f85149'; ctx.shadowBlur=5;
      ctx.beginPath(); ctx.moveTo(px,py); ctx.lineTo(tipX,tipY); ctx.stroke();
      const ah=7;
      ctx.beginPath(); ctx.moveTo(tipX,tipY);
      ctx.lineTo(tipX-ah*Math.cos(dirAng-0.4)*sign, tipY-ah*Math.sin(dirAng-0.4)*sign);
      ctx.lineTo(tipX-ah*Math.cos(dirAng+0.4)*sign, tipY-ah*Math.sin(dirAng+0.4)*sign);
      ctx.closePath(); ctx.fill();
      ctx.font='bold 11px "Consolas",monospace'; ctx.fillStyle='#ff7b72'; ctx.textAlign='center';
      const labelX=tipX+(ld.direction==='x'? (sign>0?-5:5) : 12);
      const labelY=tipY+(ld.direction==='y'? (sign>0?8:-8) : -6);
      ctx.fillText(ld.value+' kN', labelX, labelY);
      ctx.restore();
    } else if(ld.type==='moment'){
      const t=ld.pos/mem.L;
      const px=sa.x+(sb.x-sa.x)*t, py=sa.y+(sb.y-sa.y)*t;
      ctx.save();
      ctx.strokeStyle='#a371f7'; ctx.fillStyle='#a371f7'; ctx.lineWidth=2;
      ctx.shadowColor='#a371f7'; ctx.shadowBlur=5;
      const r=16, sign=ld.value>=0?1:-1;
      ctx.beginPath(); ctx.arc(px,py,r,-Math.PI/2,Math.PI/2, sign<0); ctx.stroke();
      const ea=sign<0?Math.PI/2:-Math.PI/2;
      const ex=px+r*Math.cos(ea), ey=py+r*Math.sin(ea);
      ctx.beginPath(); ctx.moveTo(ex,ey); ctx.lineTo(ex-sign*7,ey-sign*5); ctx.lineTo(ex+sign*7,ey-sign*5); ctx.closePath(); ctx.fill();
      ctx.font='bold 10.5px "Consolas",monospace'; ctx.fillStyle='#bc8cff'; ctx.textAlign='center';
      ctx.fillText(ld.value+' kN·m', px, py-r-6);
      ctx.restore();
    }
  });
}

function drawDiagram(mf, type, color){
  const scale = type==='M'? 25 : type==='V'? 10 : 10;
  const sign = type==='M'? -1 : 1; // M画在受拉侧(下方/右侧)
  
  MEMBERS.forEach(mem => {
    const f = mf[mem.id]; if(!f) return;
    const sa=toS(NODES[mem.nodeA].x, NODES[mem.nodeA].y);
    const sb=toS(NODES[mem.nodeB].x, NODES[mem.nodeB].y);
    const angle=Math.atan2(sb.y-sa.y, sb.x-sa.x);
    const nx=-Math.sin(angle), ny=Math.cos(angle); // 垂直方向(指向左侧)
    const data=f[type]; if(!data) return;
    const n=data.length;
    
    ctx.save();
    ctx.fillStyle=color; ctx.strokeStyle=color; ctx.lineWidth=2;
    ctx.shadowColor=color; ctx.shadowBlur=6; ctx.globalAlpha=0.15;
    
    // 填充
    ctx.beginPath();
    for(let i=0;i<n;i++){
      const t=i/(n-1);
      const bx=sa.x+(sb.x-sa.x)*t, by=sa.y+(sb.y-sa.y)*t;
      const v=data[i]*scale*sign;
      const px=bx+nx*v, py=by+ny*v;
      if(i===0) ctx.moveTo(bx,by);
      ctx.lineTo(px,py);
    }
    for(let i=n-1;i>=0;i--){
      const t=i/(n-1);
      const bx=sa.x+(sb.x-sa.x)*t, by=sa.y+(sb.y-sa.y)*t;
      ctx.lineTo(bx,by);
    }
    ctx.closePath(); ctx.fill();
    
    // 轮廓
    ctx.globalAlpha=1; ctx.lineWidth=2.5;
    ctx.beginPath();
    for(let i=0;i<n;i++){
      const t=i/(n-1);
      const bx=sa.x+(sb.x-sa.x)*t, by=sa.y+(sb.y-sa.y)*t;
      const v=data[i]*scale*sign;
      const px=bx+nx*v, py=by+ny*v;
      if(i===0) ctx.moveTo(px,py); else ctx.lineTo(px,py);
    }
    ctx.stroke();
    
    // 标注最大值
    let mi=0, mv=Math.abs(data[0]);
    for(let i=1;i<n;i++) if(Math.abs(data[i])>mv){mv=Math.abs(data[i]);mi=i;}
    const mt=mi/(n-1);
    const mbx=sa.x+(sb.x-sa.x)*mt, mby=sa.y+(sb.y-sa.y)*mt;
    const mpx=mbx+nx*data[mi]*scale*sign, mpy=mby+ny*data[mi]*scale*sign;
    ctx.font='bold 11px "Consolas",monospace'; ctx.fillStyle=color; ctx.textAlign='center';
    const off = data[mi]>=0? -16 : 20;
    ctx.fillText(type+'='+data[mi].toFixed(1), mpx, mpy+off);
    ctx.restore();
  });
}

function drawReactions(reactions){
  ctx.save();
  NODE_IDS.forEach(id => {
    const r=reactions[id]; if(!r) return;
    const s=toS(NODES[id].x, NODES[id].y);
    
    if(Math.abs(r.fx)>0.01){
      const dir=r.fx>=0?1:-1;
      ctx.strokeStyle='#3fb950'; ctx.fillStyle='#3fb950'; ctx.lineWidth=2.5;
      ctx.beginPath(); ctx.moveTo(s.x,s.y); ctx.lineTo(s.x+dir*32,s.y); ctx.stroke();
      ctx.beginPath(); ctx.moveTo(s.x+dir*32,s.y); ctx.lineTo(s.x+dir*25,s.y-5); ctx.lineTo(s.x+dir*25,s.y+5); ctx.closePath(); ctx.fill();
      ctx.font='bold 10.5px "Consolas",monospace'; ctx.fillStyle='#3fb950'; ctx.textAlign='left';
      ctx.fillText('Fx='+Math.abs(r.fx).toFixed(2)+'kN', s.x+dir*36, s.y-8);
    }
    if(Math.abs(r.fy)>0.01){
      const dir=r.fy>=0?1:-1;
      ctx.strokeStyle='#58a6ff'; ctx.fillStyle='#58a6ff'; ctx.lineWidth=2.5;
      ctx.beginPath(); ctx.moveTo(s.x,s.y); ctx.lineTo(s.x,s.y+dir*32); ctx.stroke();
      ctx.beginPath(); ctx.moveTo(s.x,s.y+dir*32); ctx.lineTo(s.x-5,s.y+dir*25); ctx.lineTo(s.x+5,s.y+dir*25); ctx.closePath(); ctx.fill();
      ctx.font='bold 10.5px "Consolas",monospace'; ctx.fillStyle='#58a6ff'; ctx.textAlign='left';
      ctx.fillText('Fy='+Math.abs(r.fy).toFixed(2)+'kN', s.x+8, s.y+dir*36);
    }
  });
  ctx.restore();
}

// ========== UI ==========
function updateLoadsList(){
  const list=document.getElementById('loads-list'); list.innerHTML='';
  loads.forEach(ld => {
    const mem=MEMBERS.find(m=>m.id===ld.member);
    const tl=ld.type==='udl'?'均布荷载':ld.type==='point'?'集中力':'集中力偶';
    const dl=ld.direction==='x'?'水平':ld.direction==='y'?'竖向':'';
    const div=document.createElement('div'); div.className='load-item';
    div.innerHTML=`<div class="load-hdr"><span class="load-type">${tl}</span><button class="del-btn" data-id="${ld.id}">✕</button></div>
      <label>${ld.member}杆 (${mem.L}m) · ${ld.type==='udl'?ld.pos1+'→'+ld.pos2+'m':'位置:'+ld.pos+'m'} ${dl}</label>
      <label>值: ${ld.value} ${ld.type==='udl'?'kN/m':ld.type==='moment'?'kN·m':'kN'}</label>`;
    div.querySelector('.del-btn').onclick=()=>{
      loads=loads.filter(l=>l.id!==ld.id); solved=false; solution=null;
      updateLoadsList(); render(); updateStatus('已删除荷载');
    };
    list.appendChild(div);
  });
}

function updateMemberSelect(){
  const sel=document.getElementById('new-load-pos'); sel.innerHTML='';
  MEMBERS.forEach(mem => {
    const opt=document.createElement('option'); opt.value=mem.id;
    const cnt=loads.filter(l=>l.member===mem.id).length;
    opt.textContent=`${mem.id}杆 (${mem.L}m)${cnt?` · ${cnt}个荷载`:''}`;
    sel.appendChild(opt);
  });
  updateLoadUI();
}

function updateLoadUI(){
  const type=document.getElementById('new-load-type').value;
  const dirRow=document.getElementById('dir-row');
  const posRow=document.getElementById('pos-row');
  if(type==='point'){ dirRow.style.display='flex'; posRow.style.display='flex'; }
  else if(type==='moment'){ dirRow.style.display='none'; posRow.style.display='flex'; }
  else { dirRow.style.display='none'; posRow.style.display='none'; }
}

function updateResults(){
  if(!solution) return;
  // 反力
  const rx=document.getElementById('reactions'); rx.innerHTML='';
  NODE_IDS.forEach(id => {
    const r=solution.reactions[id];
    if(Math.abs(r.fx)<0.01 && Math.abs(r.fy)<0.01) return;
    const d=document.createElement('div'); d.className='result-block';
    let h=`<div class="rb-title">节点 ${id}</div>`;
    if(Math.abs(r.fx)>0.01) h+=`<div><span class="rb-val" style="color:#3fb950">Fx = ${Math.abs(r.fx).toFixed(2)}</span> <span class="rb-unit">kN ${r.fx>=0?'→':'←'}</span></div>`;
    if(Math.abs(r.fy)>0.01) h+=`<div><span class="rb-val" style="color:#58a6ff">Fy = ${Math.abs(r.fy).toFixed(2)}</span> <span class="rb-unit">kN ${r.fy>=0?'↑':'↓'}</span></div>`;
    d.innerHTML=h; rx.appendChild(d);
  });
  // 杆件内力
  const mr=document.getElementById('member-results'); mr.innerHTML='';
  MEMBERS.forEach(mem => {
    const f=solution.memberForces[mem.id]; if(!f) return;
    const mM=Math.max(...f.M.map(Math.abs));
    const mV=Math.max(...f.V.map(Math.abs));
    const mN=Math.max(...f.N.map(Math.abs));
    const d=document.createElement('div'); d.className='member-result';
    d.innerHTML=`<div class="mr-hdr"><span class="mr-name">${mem.id}杆 (L=${mem.L}m)</span></div>
      <div class="mr-row"><span>M_max: <span class="v">${mM.toFixed(2)} kN·m</span></span>
      <span>V_max: <span class="v">${mV.toFixed(2)} kN</span></span>
      <span>N_max: <span class="v">${mN.toFixed(2)} kN</span></span></div>
      <div class="mr-row" style="margin-top:4px;color:#6e7681;font-size:10px;">
      ①: M=${f.start.M.toFixed(1)}, V=${f.start.V.toFixed(1)}, N=${f.start.N.toFixed(1)}
      ②: M=${f.end.M.toFixed(1)}, V=${f.end.V.toFixed(1)}, N=${f.end.N.toFixed(1)}</div>`;
    mr.appendChild(d);
  });
}

function updateStatus(msg){ document.getElementById('status').innerHTML=msg; }

// ========== 事件 ==========
document.getElementById('btn-solve').onclick=() => {
  try {
    solution=solve(); solved=true;
    updateResults(); render();
    updateStatus('求解完成 ✓ 切到「计算结果」查看详情');
  } catch(e) {
    updateStatus('<span style="color:#f85149">求解失败: '+e.message+'</span>');
  }
};
document.getElementById('btn-show-m').onclick=e=>{showM=!showM;e.target.classList.toggle('active',showM);render();};
document.getElementById('btn-show-v').onclick=e=>{showV=!showV;e.target.classList.toggle('active',showV);render();};
document.getElementById('btn-show-n').onclick=e=>{showN=!showN;e.target.classList.toggle('active',showN);render();};
document.getElementById('btn-show-frames').onclick=e=>{showFrame=!showFrame;e.target.classList.toggle('active',showFrame);render();};
document.getElementById('btn-reset').onclick=()=>{
  loads=[
    {id:1,type:'udl',member:'CD',value:10,pos1:0,pos2:6},
    {id:2,type:'point',member:'DE',pos:2,value:9,direction:'y'},
    {id:3,type:'point',member:'BD',pos:2,value:6,direction:'x'}
  ];
  solved=false; solution=null;
  updateLoadsList(); render();
  updateStatus('已恢复默认荷载');
};
document.querySelectorAll('.tab').forEach(t => {
  t.onclick=() => {
    document.querySelectorAll('.tab').forEach(x=>x.classList.remove('active'));
    t.classList.add('active');
    const tab=t.dataset.tab;
    document.getElementById('tab-loads').style.display=tab==='loads'?'block':'none';
    document.getElementById('tab-results').style.display=tab==='results'?'block':'none';
    if(tab==='results') updateResults();
  };
});
document.getElementById('new-load-type').onchange=updateLoadUI;
document.getElementById('btn-add-load').onclick=() => {
  const type=document.getElementById('new-load-type').value;
  const mem=document.getElementById('new-load-pos').value;
  const val=parseFloat(document.getElementById('new-load-value').value);
  if(isNaN(val)){ updateStatus('<span style="color:#f85149">请输入有效数值</span>'); return; }
  const newId=loads.length?Math.max(...loads.map(l=>l.id))+1:1;
  if(type==='udl'){
    const m=MEMBERS.find(x=>x.id===mem);
    loads.push({id:newId,type:'udl',member:mem,value:val,pos1:0,pos2:m.L});
  } else if(type==='point'){
    const dir=document.getElementById('new-load-dir').value;
    const pos=parseFloat(document.getElementById('new-load-pos-val').value)||1;
    loads.push({id:newId,type:'point',member:mem,pos,value:val,direction:dir});
  } else {
    const pos=parseFloat(document.getElementById('new-load-pos-val').value)||1;
    loads.push({id:newId,type:'moment',member:mem,pos,value:val});
  }
  solved=false; solution=null;
  updateLoadsList(); render();
  updateStatus('已添加荷载 ✓ 点击「求解」');
};

// 初始化
window.addEventListener('load', () => { setTimeout(resize, 50); });
resize(); updateLoadsList(); updateMemberSelect(); render();
