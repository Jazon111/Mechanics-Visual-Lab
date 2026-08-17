// 测试分层求解：附属/基本结构识别 + 先附属后基本 + 截面法
const E = 2e8, A_sec = 0.01, I_sec = 0.0001;
let NODE_IDS = [], NODES = {}, MEMBERS = [], SUPPORT_TYPES = {}, hingeNodes = {};
let DOF_TOTAL = 0, dofMap = {};
let loads = [];

function buildFromLengths(lengths) {
  NODE_IDS = []; NODES = {}; MEMBERS = [];
  const base = 'A'.charCodeAt(0);
  for(let i=0;i<lengths.length+1 && i<14;i++) NODE_IDS.push(String.fromCharCode(base+i));
  let x = 0;
  NODE_IDS.forEach((id,i) => {
    NODES[id] = {id, x, y:0};
    if(i < lengths.length){
      MEMBERS.push({id:NODE_IDS[i]+NODE_IDS[i+1], nodeA:id, nodeB:NODE_IDS[i+1], L:lengths[i], angle:0});
      x += lengths[i];
    }
  });
  DOF_TOTAL = NODE_IDS.length * 3;
  dofMap = {};
  NODE_IDS.forEach((id, i) => { dofMap[id] = [i*3, i*3+1, i*3+2]; });
  SUPPORT_TYPES = {};
  NODE_IDS.forEach(id => { SUPPORT_TYPES[id] = 'none'; });
  hingeNodes = {};
}
function getSupportConstraints(t) {
  return t==='fixed'?[true,true,true]:t==='pinned'?[true,true,false]:t==='roller'?[false,true,false]:t==='slider'?[true,false,true]:[false,false,false];
}
function getGlobalForce(ld){const P=ld.value;return ld.direction==='down'?{Fx:0,Fy:-P}:{Fx:0,Fy:0};}
function getMomentSign(ld){return ld.direction==='cw'?-1:1;}
function gaussSolve(A,b){const n=A.length,M=A.map((r,i)=>[...r,b[i]]);for(let i=0;i<n;i++){let mr=i;for(let k=i+1;k<n;k++)if(Math.abs(M[k][i])>Math.abs(M[mr][i]))mr=k;[M[i],M[mr]]=[M[mr],M[i]];const p=M[i][i];if(Math.abs(p)<1e-14)continue;for(let k=i;k<=n;k++)M[i][k]/=p;for(let k=0;k<n;k++){if(k===i)continue;const f=M[k][i];for(let j=i;j<=n;j++)M[k][j]-=f*M[i][j];}}return M.map(r=>r[n]);}

// === 核心函数（从HTML复制） ===
function identifyLayers() {
  const layers = [];
  let startIdx = 0;
  for(let i=1; i<NODE_IDS.length; i++) {
    const id = NODE_IDS[i];
    if(hingeNodes[id] || i === NODE_IDS.length-1) {
      const nodes = NODE_IDS.slice(startIdx, i+1);
      const members = MEMBERS.filter(m => nodes.includes(m.nodeA) && nodes.includes(m.nodeB));
      const vertSupports = nodes.filter(n => getSupportConstraints(SUPPORT_TYPES[n] || 'none')[1]);
      const isBasic = vertSupports.length >= 2;
      layers.push({id:`L${layers.length+1}`,startNode:NODE_IDS[startIdx],endNode:id,nodes,members,vertSupports,isBasic,type:isBasic?'基本结构':'附属结构'});
      startIdx = i;
    }
  }
  return layers;
}
function collectLayerForces(layer, hingeTransfer) {
  const ext = [];
  layer.members.forEach(mem => {
    const gxA = NODES[mem.nodeA].x;
    loads.filter(ld => ld.member === mem.id).forEach(ld => {
      if(ld.type==='udl'){const w=ld.value,a=ld.pos1,b=ld.pos2;if(b>a)ext.push({Fx:0,Fy:-w*(b-a),M:0,gx:gxA+(a+b)/2,src:'udl'});}
      else if(ld.type==='point'){const {Fx,Fy}=getGlobalForce(ld);ext.push({Fx,Fy,M:0,gx:gxA+ld.pos,src:'point'});}
      else if(ld.type==='moment'){ext.push({Fx:0,Fy:0,M:ld.value*getMomentSign(ld),gx:gxA+ld.pos,src:'moment'});}
    });
  });
  layer.nodes.forEach(n => {
    if(hingeNodes[n] && hingeTransfer[n] !== undefined)
      ext.push({Fx:0,Fy:hingeTransfer[n],M:0,gx:NODES[n].x,src:'hinge-transfer'});
  });
  return ext;
}
function solveLayerEquilibrium(layer, reactions, hingeTransfer) {
  const ext = collectLayerForces(layer, hingeTransfer);
  const unk = [];
  layer.nodes.forEach(n => {
    if(!layer.isBasic && (n === layer.startNode || n === layer.endNode) && hingeNodes[n]) return;
    if(hingeNodes[n] && hingeTransfer[n] !== undefined) {
      const sc0 = getSupportConstraints(SUPPORT_TYPES[n] || 'none');
      if(!sc0[0] && !sc0[1] && !sc0[2]) return;
    }
    const sc = getSupportConstraints(SUPPORT_TYPES[n] || 'none');
    if(sc[0]) unk.push({node:n,comp:'fx'});
    if(sc[1]) unk.push({node:n,comp:'fy'});
    if(sc[2]) unk.push({node:n,comp:'m'});
  });
  if(!layer.isBasic) {
    [layer.startNode, layer.endNode].forEach(n => {
      if(hingeNodes[n] && hingeTransfer[n] === undefined)
        unk.push({node:n,comp:'fy',isHinge:true});
    });
  }
  const n = unk.length;
  if(n === 0) return {ext, unk, values:[]};
  const A = Array.from({length:n},()=>Array(n).fill(0));
  const b = Array(n).fill(0);
  const x0 = NODES[layer.startNode].x, xR = NODES[layer.endNode].x;
  let eq = 0;
  if(unk.some(u=>u.comp==='fx')) {
    unk.forEach((u,i)=>{if(u.comp==='fx')A[eq][i]=1;});
    b[eq] = -ext.reduce((s,f)=>s+f.Fx,0);
    eq++;
  }
  if(eq < n) {
    unk.forEach((u,i)=>{if(u.comp==='fy')A[eq][i]=1;});
    b[eq] = -ext.reduce((s,f)=>s+f.Fy,0);
    eq++;
  }
  if(eq < n) {
    unk.forEach((u,i)=>{const xu=NODES[u.node].x;if(u.comp==='fy')A[eq][i]=(x0-xu);if(u.comp==='m')A[eq][i]=-1;});
    b[eq] = -ext.reduce((s,f)=>s+f.Fy*(x0-f.gx)-f.M,0);
    eq++;
  }
  layer.nodes.forEach(nd => {
    if(eq >= n) return;
    if(hingeNodes[nd] && nd !== layer.startNode && nd !== layer.endNode) {
      const xh = NODES[nd].x;
      unk.forEach((u,i)=>{const xu=NODES[u.node].x;if(xu<xh+1e-9){if(u.comp==='fy')A[eq][i]=(xh-xu);if(u.comp==='m')A[eq][i]=-1;}});
      b[eq] = -ext.filter(f=>f.gx<xh+1e-9).reduce((s,f)=>s+f.Fy*(xh-f.gx)-f.M,0);
      eq++;
    }
  });
  for(; eq < n; eq++) {
    unk.forEach((u,i)=>{const xu=NODES[u.node].x;if(u.comp==='fy')A[eq][i]=(xR-xu);if(u.comp==='m')A[eq][i]=-1;});
    b[eq] = -ext.reduce((s,f)=>s+f.Fy*(xR-f.gx)-f.M,0);
  }
  const x = gaussSolve(A, b);
  const values = unk.map((u,i) => isFinite(x[i]) ? x[i] : 0);
  unk.forEach((u,i) => {
    const val = values[i];
    if(u.isHinge) { hingeTransfer[u.node] = -val; }
    else {
      if(u.comp==='fx') reactions[u.node].fx += val;
      if(u.comp==='fy') reactions[u.node].fy += val;
      if(u.comp==='m') reactions[u.node].m += val;
    }
  });
  return {ext, unk, values};
}
function solveStaticallyDeterminate() {
  const layers = identifyLayers();
  const reactions = {};
  NODE_IDS.forEach(id => { reactions[id] = {fx:0,fy:0,m:0}; });
  const hingeTransfer = {};
  const order = [...layers.filter(l=>!l.isBasic), ...layers.filter(l=>l.isBasic)];
  order.forEach(layer => { solveLayerEquilibrium(layer, reactions, hingeTransfer); });
  return {reactions, layers, order, hingeTransfer};
}

// === 截面法（简化版，只测关键点） ===
function computeM_at(gxSec, reactions) {
  let sumM = 0;
  NODE_IDS.forEach(id => {
    const r = reactions[id];
    const gx = NODES[id].x;
    if(gx < gxSec - 1e-9) sumM += r.fy * (gxSec - gx);
    else if(Math.abs(gx - gxSec) < 1e-9 && gx === NODES[NODE_IDS[0]].x) sumM += r.fy * (gxSec - gx);
  });
  loads.forEach(ld => {
    const mem = MEMBERS.find(m=>m.id===ld.member); if(!mem) return;
    const gxA = NODES[mem.nodeA].x;
    if(ld.type==='udl'){const w=ld.value,p1=ld.pos1,p2=ld.pos2;const gx1=gxA+p1,gx2=gxA+p2;if(gxSec>gx1){const effEnd=Math.min(gxSec,gx2);const len=effEnd-gx1;if(len>0){const Fy=-w*len;const pos=gx1+len/2;sumM+=Fy*(gxSec-pos);}}}
    else if(ld.type==='point'){const gx=gxA+ld.pos;if(gx<gxSec-1e-9){const {Fy}=getGlobalForce(ld);sumM+=Fy*(gxSec-gx);}}
    else if(ld.type==='moment'){const gx=gxA+ld.pos;if(gx<gxSec-1e-9)sumM-=ld.value*getMomentSign(ld);}
  });
  return sumM;
}

let passCount = 0, failCount = 0;
function assert(name, got, expected, tol=0.05) {
  const diff = Math.abs(got - expected);
  const ok = diff < tol;
  if(ok) { console.log(`  ✓ ${name}: got=${got.toFixed(3)}, expected=${expected.toFixed(3)}`); passCount++; }
  else { console.log(`  ✗ ${name}: got=${got.toFixed(3)}, expected=${expected.toFixed(3)} FAIL`); failCount++; }
}

// ===== 例3: 简支梁L=6 UDL w=10 =====
console.log('\n==== 例3: 简支梁 L=6m UDL w=10 ====');
buildFromLengths([6]);
SUPPORT_TYPES['A']='pinned'; SUPPORT_TYPES['B']='roller';
loads=[{id:1,type:'udl',member:'AB',value:10,pos1:0,pos2:6}];
const sol3 = solveStaticallyDeterminate();
console.log('  Layers:', sol3.layers.map(l=>`${l.id}[${l.startNode}-${l.endNode}]${l.type}`));
console.log('  Order:', sol3.order.map(l=>l.id));
console.log('  Reactions: A=', sol3.reactions.A, 'B=', sol3.reactions.B);
assert('R_Ay=30', sol3.reactions.A.fy, 30);
assert('R_By=30', sol3.reactions.B.fy, 30);
assert('M@3m=45', computeM_at(3, sol3.reactions), 45);
assert('M@0=0', computeM_at(0, sol3.reactions), 0);
assert('M@6=0', computeM_at(6, sol3.reactions), 0);

// ===== 例2: 简支6m+悬臂2m =====
console.log('\n==== 例2: 简支6m+悬臂2m ====');
buildFromLengths([6, 2]);
SUPPORT_TYPES['A']='pinned'; SUPPORT_TYPES['B']='roller';
loads=[
  {id:1,type:'udl',member:'AB',value:10,pos1:0,pos2:6},
  {id:2,type:'point',member:'BC',pos:2,value:30,direction:'down'}
];
const sol2 = solveStaticallyDeterminate();
console.log('  Layers:', sol2.layers.map(l=>`${l.id}[${l.startNode}-${l.endNode}]${l.type}`));
console.log('  Reactions: A=', sol2.reactions.A, 'B=', sol2.reactions.B);
assert('R_Ay=20', sol2.reactions.A.fy, 20);
assert('R_By=70', sol2.reactions.B.fy, 70);
assert('M@6m(B点)=-60', computeM_at(6, sol2.reactions), -60);
assert('M@2m=+20', computeM_at(2, sol2.reactions), 20);
assert('M@8m(C点)=0', computeM_at(8, sol2.reactions), 0);

// ===== 例1: 多跨静定梁（修改后F=roller） =====
console.log('\n==== 例1: 多跨静定梁 A-B-C(铰)-D(铰+可动)-E-F(可动) ====');
buildFromLengths([4, 1, 3, 1, 4]);
SUPPORT_TYPES['A']='pinned'; SUPPORT_TYPES['B']='roller'; SUPPORT_TYPES['D']='roller'; SUPPORT_TYPES['F']='roller';
hingeNodes['C']=true; hingeNodes['D']=true;
loads = [
  {id:1,type:'point',member:'AB',pos:1,value:8,direction:'down'},
  {id:2,type:'point',member:'AB',pos:3,value:8,direction:'down'},
  {id:3,type:'moment',member:'CD',pos:1,value:24,direction:'ccw'},
  {id:4,type:'udl',member:'EF',value:4,pos1:0,pos2:4}
];
const sol1 = solveStaticallyDeterminate();
console.log('  Layers:');
sol1.layers.forEach(l => console.log(`    ${l.id} [${l.startNode}→${l.endNode}] ${l.type} (${l.vertSupports.length}竖向支座)`));
console.log('  Solve order:', sol1.order.map(l=>`${l.id}(${l.type})`));
console.log('  Hinge transfers:', sol1.hingeTransfer);
console.log('  Reactions:');
NODE_IDS.forEach(id => { const r=sol1.reactions[id]; if(SUPPORT_TYPES[id]!=='none') console.log(`    ${id}: fy=${r.fy.toFixed(3)}`); });

// 验证整体平衡
let sumFy = 0;
NODE_IDS.forEach(id => { sumFy += sol1.reactions[id].fy; });
sumFy += (-8) + (-8) + (-4*4); // 外荷载
assert('ΣFy=0', sumFy, 0, 0.01);

// 验证铰C弯矩=0
const MC = computeM_at(NODES['C'].x, sol1.reactions);
assert('M_C(铰)=0', MC, 0, 0.1);

// 验证铰D弯矩=0
const MD = computeM_at(NODES['D'].x, sol1.reactions);
assert('M_D(铰)=0', MD, 0, 0.1);

// 验证关键截面弯矩
// B点(x=4): M = R_Ay*4 - P1*3 - P2*1 = 6*4 - 8*3 - 8*1 = -8
assert('M_B=-8', computeM_at(4, sol1.reactions), -8, 0.1);

// 集中力偶左侧(x=5.9): M = 6*5.9+18*1.9-8*4.9-8*2.9 = 35.4+34.2-39.2-23.2 = 7.2
// (D在x=8>5.9, 不计入左侧; 力偶在x=6>5.9, 也不计入)
assert('M_左 Couple=7.2', computeM_at(5.9, sol1.reactions), 7.2, 0.1);

// 集中力偶右侧(x=6.1): M = 7.2 + V*0.2 - M_couple = 7.2 + 8*0.2 - 24 = -15.2
// (V = R_A+R_B-P1-P2 = 6+18-8-8 = 8; 力偶CCW使弯矩图向下跳跃24)
assert('M_右 Couple=-15.2', computeM_at(6.1, sol1.reactions), -15.2, 0.1);

// 验证力偶跳跃: M_right - M_left = -24 (CCW力偶使弯矩向下跳跃)
const M_jump = computeM_at(6.1, sol1.reactions) - computeM_at(5.9, sol1.reactions) - 8*0.2;
assert('Couple jump=-24', M_jump, -24, 0.1);

// E点(x=9, UDL起点): M = 6*9+18*5-1.6*1-8*8-8*6-24 = 54+90-1.6-64-48-24 = 6.4
// 从右侧验证: -[R_Fy*(9-13)+UDL*(9-11)] = -[9.6*(-4)+(-16)*(-2)] = -[-38.4+32] = 6.4
assert('M_E@x=9', computeM_at(9, sol1.reactions), 6.4, 0.5);

console.log('\n==== Summary: PASS=' + passCount + ' FAIL=' + failCount + ' ====');
process.exit(failCount>0 ? 1 : 0);
