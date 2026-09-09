/*
  version: 1.0.0
  
  UPDATES(m/d/y): 
	9/6/2026 - added terminal inputs
	9/6/2026 - performance fix, replaced object vectors with inline math inside loops
	9/6/2026 - added angle forces
	9/6/2026 - fixed issue with temperatures starting at over 1000k
	9/7/2026 - fixed memory issues
	9/7/2026 - added multi-threading and amd gpu support (havent tested nvidia)
	9/8/2026 - spent the entire day fixing issues that were holding back realism and implementing new formulas, bigger update here for realism
	9/8/2026 - CIC particle-mesh Ewald reciprocal solve plus Ewald real space.  Triclinic cells not supported.
	9/8/2026 - posting code online
	
  TO DO:
	- add comments and format code to make it readable
	- be able to specify trajectory output location and name
	
	later:
		- add xml and other files to get more data about protein
*/

const fs = require('fs');
const readline = require('readline');
const os = require('os');

// system units
// CHARMM standard
const config = {
  // Constrained X-H bonds to permit a 1 fs step
  dt: 1.0, // time in femtoseconds
  thermostatRelaxationRate: 100.0, // relaxation rate tau (fs)
  gamma: 0.01, // friction (1/fs) as 1/thermostatRelaxationRate
  targetTemp: 300, // kelvin
  kB: 0.0019872041, // kcal/(mol*k)
  fta: 4.184e-4, // kcal/(mol*A*amu) to A/fs^2
  coulomb: 332.06371, // kcal*A/(mol*e^2)
  cutoff: 10.0,
  switchDistance: 8.0,
  neighborSkin: 2.0,
  oneFourLJScale: 0.5,
  oneFourCoulombScale: 1 / 1.2,
  constraintTolerance: 1e-6,
  constraintIterations: 100,
  pme: true,
  pmeGridSpacing: 1.0, // Å
  ewaldAlpha: 0.35 // Å^-1 erfc(alpha * cutoff)
};

const masses = {
  H: 1.008, He: 4.0026, Li: 6.94, Be: 9.0122, B: 10.81, C: 12.011, N: 14.007, O: 15.999, F: 18.998, Ne: 20.180,
  Na: 22.990, Mg: 24.305, Al: 26.982, Si: 28.085, P: 30.974, S: 32.06, Cl: 35.45, Ar: 39.95, K: 39.098, Ca: 40.078,
  Sc: 44.956, Ti: 47.867, V: 50.942, Cr: 51.996, Mn: 54.938, Fe: 55.845, Co: 58.933, Ni: 58.693, Cu: 63.546, Zn: 65.38,
  Ga: 69.723, Ge: 72.630, As: 74.922, Se: 78.971, Br: 79.904, Kr: 83.798, Rb: 85.468, Sr: 87.62, Y: 88.906, Zr: 91.224,
  Nb: 92.906, Mo: 95.95, Tc: 98, Ru: 101.07, Rh: 102.91, Pd: 106.42, Ag: 107.87, Cd: 112.41, In: 114.82, Sn: 118.71,
  Sb: 121.76, Te: 127.60, I: 126.90, Xe: 131.29, Cs: 132.91, Ba: 137.33, La: 138.91, Ce: 140.12, Pr: 140.91, Nd: 144.24,
  Pm: 145, Sm: 150.36, Eu: 151.96, Gd: 157.25, Tb: 158.93, Dy: 162.50, Ho: 164.93, Er: 167.26, Tm: 168.93, Yb: 173.05,
  Lu: 174.97, Hf: 178.49, Ta: 180.95, W: 183.84, Re: 186.21, Os: 190.23, Ir: 192.22, Pt: 195.08, Au: 196.97, Hg: 200.59,
  Tl: 204.38, Pb: 207.2, Bi: 208.98, Po: 209, At: 210, Rn: 222, Fr: 223, Ra: 226, Ac: 227, Th: 232.04,
  Pa: 231.04, U: 238.03, Np: 237, Pu: 244, Am: 243, Cm: 247, Bk: 247, Cf: 251, Es: 252, Fm: 257,
  Md: 258, No: 259, Lr: 262, Rf: 267, Db: 270, Sg: 271, Bh: 270, Hs: 277, Mt: 276, Ds: 281,
  Rg: 282, Cn: 285, Nh: 286, Fl: 289, Mc: 290, Lv: 293, Ts: 294, Og: 294,
  undefDef: 12.0
};

const ff = {
  lj: {
    H:  {eps:0.0300,sig:2.50 },
    C:  {eps:0.0660,sig:3.50 },
    N:  {eps:0.1700,sig:3.25 },
    O:  {eps:0.1521,sig:3.15061 },
    S:  {eps:0.2500,sig:3.60 },
    P:  {eps:0.2000,sig:3.74 },
    F:  {eps:0.0610,sig:3.12 },
    Cl: {eps:0.2270,sig:3.47 },
    Br: {eps:0.3200,sig:3.68 },
    I:  {eps:0.4000,sig:4.00 },
    Na: {eps:0.0028,sig:3.33 },
    K:  {eps:0.0003,sig:4.73 },
    Mg: {eps:0.8750,sig:1.41 },
    Ca: {eps:0.1200,sig:2.87 },
    undefDef: {eps:0.1,sig:3.0 }
  },
  bonds: {
    'C-H': { k: 340.0, r0: 1.090 },
    'H-O': { k: 450.0, r0: 0.9572 },
    'C-C': { k: 310.0, r0: 1.526 },
    'C-O': { k: 320.0, r0: 1.410 },
    'H-N': { k: 434.0, r0: 1.010 },
    'C-N': { k: 320.0, r0: 1.470 },
    'C-S': { k: 250.0, r0: 1.810 },
    'S-S': { k: 166.0, r0: 2.040 },
    'S-H': { k: 274.0, r0: 1.336 },
    'C-F': { k: 367.0, r0: 1.350 },
    'C-Cl':{ k: 300.0, r0: 1.770 },
    'C-Br':{ k: 280.0, r0: 1.940 },
    'C-I': { k: 230.0, r0: 2.140 },
    'P-O': { k: 400.0, r0: 1.610 },
    undefDef: { k: 300, r0: 1.5 }
  },
  angles: {
    'H-C-H': { k: 39.5, t0: 109.5 },
    'H-O-H': { k: 55.0, t0: 104.52 },
    'H-C-C': { k: 50.0, t0: 109.5 },
    'C-C-C': { k: 40.0, t0: 109.5 },
    'C-C-O': { k: 50.0, t0: 109.5 },
    'C-C-N': { k: 50.0, t0: 109.5 },
    'C-N-C': { k: 50.0, t0: 109.5 },
    'C-O-C': { k: 60.0, t0: 109.5 },
    'O-P-O': { k: 100.0, t0: 109.5 },
    'C-S-C': { k: 62.0, t0: 98.9 },
    undefDef: { k: 50, t0: 109.5 }
  },
  dihedrals: {
    undefDef: { k: 1.0, n: 1, delta: 0.0 }
  }
};

const backboneBondPairs = [['N','CA'],['CA','C'],['C','O'],['C','OXT']];

const residueBondPairs = {
  ALA: [['CA','CB']],
  ARG: [['CA','CB'],['CB','CG'],['CG','CD'],['CD','NE'],['NE','CZ'],['CZ','NH1'],['CZ','NH2']],
  ASN: [['CA','CB'],['CB','CG'],['CG','OD1'],['CG','ND2']],
  ASP: [['CA','CB'],['CB','CG'],['CG','OD1'],['CG','OD2']],
  CYS: [['CA','CB'],['CB','SG']],
  GLN: [['CA','CB'],['CB','CG'],['CG','CD'],['CD','OE1'],['CD','NE2']],
  GLU: [['CA','CB'],['CB','CG'],['CG','CD'],['CD','OE1'],['CD','OE2']],
  GLY: [],
  HIS: [['CA','CB'],['CB','CG'],['CG','ND1'],['ND1','CE1'],['CE1','NE2'],['NE2','CD2'],['CD2','CG']],
  ILE: [['CA','CB'],['CB','CG1'],['CB','CG2'],['CG1','CD1']],
  LEU: [['CA','CB'],['CB','CG'],['CG','CD1'],['CG','CD2']],
  LYS: [['CA','CB'],['CB','CG'],['CG','CD'],['CD','CE'],['CE','NZ']],
  MET: [['CA','CB'],['CB','CG'],['CG','SD'],['SD','CE']],
  PHE: [['CA','CB'],['CB','CG'],['CG','CD1'],['CD1','CE1'],['CE1','CZ'],['CZ','CE2'],['CE2','CD2'],['CD2','CG']],
  PRO: [['CA','CB'],['CB','CG'],['CG','CD'],['CD','N']],
  SER: [['CA','CB'],['CB','OG']],
  THR: [['CA','CB'],['CB','OG1'],['CB','CG2']],
  TRP: [['CA','CB'],['CB','CG'],['CG','CD1'],['CD1','NE1'],['NE1','CE2'],['CE2','CD2'],['CD2','CG'],['CE2','CZ2'],['CZ2','CH2'],['CH2','CZ3'],['CZ3','CE3'],['CE3','CD2']],
  TYR: [['CA','CB'],['CB','CG'],['CG','CD1'],['CD1','CE1'],['CE1','CZ'],['CZ','CE2'],['CE2','CD2'],['CD2','CG'],['CZ','OH']],
  VAL: [['CA','CB'],['CB','CG1'],['CB','CG2']]
};

class Vec3 {
  constructor(x = 0, y = 0, z = 0) {
    this.x = x; this.y = y; this.z = z;
  }
  add(v) { this.x += v.x; this.y += v.y; this.z += v.z; return this; }
  sub(v) { return new Vec3(this.x-v.x, this.y-v.y, this.z-v.z); }
  scale(s) { this.x *= s; this.y *= s; this.z *= s; return this; }
  magSq() { return this.x*this.x+this.y*this.y+this.z*this.z; }
  mag() { return Math.sqrt(this.magSq()); }
  dot(v) { return this.x*v.x+this.y*v.y+this.z*v.z; }
  copy() { return new Vec3(this.x, this.y, this.z); }
}

class Atom {
  constructor(element, pos, q = 0, serial = 0, type = null) {
    this.element = element;
    this.mass = masses[element] || masses.undefDef;
    this.pos = pos;
    this.vel = new Vec3();
    this.q = q;
    this.serial = serial;
    this.type = type || element;
  }
}

class System {
  constructor(atoms, bonds, box = null) {
    this.atoms = atoms;
    this.box = box;
    this.bonds = bonds.map(([i, j]) => this._createBond(i, j));
    this.angles = this._ba(bonds);
    this.dihedrals = this._buildDihedrals(bonds);
    this.pairClasses = this._buildPairClasses();
    this.forces = atoms.map(() => new Vec3());
    this.neighborList = new NeighborList(this);
    this.pme = config.pme && this.box ? new PMESolver(this.box) : null;
  }
  _getPairKey(a, b) {
    return [a, b].sort().join('-');
  }
  _createBond(i, j) {
    const typePair = this._getPairKey(this.atoms[i].type, this.atoms[j].type);
    const elementPair = this._getPairKey(this.atoms[i].element, this.atoms[j].element);
    const params = ff.bonds[typePair] || ff.bonds[elementPair] || ff.bonds.undefDef;
    return { i, j, k: params.k, r0: params.r0 };
  }
  _ba(bonds) {
    const adj = Array.from({ length: this.atoms.length }, () => []);
    for (const [i, j] of bonds) {
      adj[i].push(j);
      adj[j].push(i);
    }
    const angles = [];
    adj.forEach((neighbors, j) => {
      for (let a = 0; a < neighbors.length; a++) {
        for (let b = a+1; b < neighbors.length; b++) {
          const i = neighbors[a];
          const k = neighbors[b];
          const typeKey = `${[this.atoms[i].type, this.atoms[k].type].sort().join('-')}-${this.atoms[j].type}`;
          const elementKey = `${[this.atoms[i].element, this.atoms[k].element].sort().join('-')}-${this.atoms[j].element}`;
          const params = ff.angles[typeKey] || ff.angles[elementKey] || ff.angles.undefDef;
          angles.push({
            i, j, k,
            kTheta: params.k,
            theta0: (params.t0*Math.PI)/180
          });
        }
      }
    });
    return angles;
  }
  _buildDihedrals(bonds) {
    const adj = Array.from({ length: this.atoms.length }, () => []);
    for (const [i, j] of bonds) {
      adj[i].push(j);
      adj[j].push(i);
    }
    const dihedrals = [];
    for (const [j, k] of bonds) {
      for (const i of adj[j]) {
        if (i === k) continue;
        for (const l of adj[k]) {
          if (l === j || l === i) continue;
          const params = ff.dihedrals.undefDef;
          dihedrals.push({
            i, j, k, l,
            k: params.k,
            n: params.n,
            delta: params.delta
          });
        }
      }
    }
    return dihedrals;
  }
  _buildPairClasses() {
    const pairClasses = new Map();
    const setClass = (i, j, value) => {
      const min = Math.min(i, j);
      const max = Math.max(i, j);
      pairClasses.set(`${min}:${max}`, value);
    };
    this.bonds.forEach(b => setClass(b.i, b.j, 1));
    this.angles.forEach(a => setClass(a.i, a.k, 1));
    this.dihedrals.forEach(d => {
      const key = `${Math.min(d.i, d.l)}:${Math.max(d.i, d.l)}`;
      if (!pairClasses.has(key)) setClass(d.i, d.l, 2);
    });
    return pairClasses;
  }
}

class NeighborList {
  constructor(sys) {
    this.sys = sys;
    this.pairs = [];
    this.reference = [];
    this.valid = false;
  }
  needsRebuild() {
    if (!this.valid) return true;
    const limit2 = (config.neighborSkin * 0.5) ** 2;
    for (let i = 0; i < this.sys.atoms.length; i++) {
      const p = this.sys.atoms[i].pos, r = this.reference[i];
      const d = pbcDistance(p.x - r.x, p.y - r.y, p.z - r.z, this.sys.box);
      if (d.dx * d.dx + d.dy * d.dy + d.dz * d.dz > limit2) return true;
    }
    return false;
  }
  rebuild() {
    const atoms = this.sys.atoms;
    const cellSize = config.cutoff + config.neighborSkin;
    const cells = new Map();
    const key = (x, y, z) => `${x},${y},${z}`;
    const dims = this.sys.box ? [
      Math.max(1, Math.ceil(this.sys.box.x / cellSize)),
      Math.max(1, Math.ceil(this.sys.box.y / cellSize)),
      Math.max(1, Math.ceil(this.sys.box.z / cellSize))
    ] : null;
    for (let i = 0; i < atoms.length; i++) {
      const p = atoms[i].pos;
      const cx = Math.floor(p.x / cellSize), cy = Math.floor(p.y / cellSize), cz = Math.floor(p.z / cellSize);
      const k = key(cx, cy, cz);
      if (!cells.has(k)) cells.set(k, []);
      cells.get(k).push(i);
    }
    const pairs = [];
    const maxR2 = cellSize * cellSize;
    for (const [cellKey, members] of cells) {
      const [cx, cy, cz] = cellKey.split(',').map(Number);
      const neighborKeys = new Set();
      for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) for (let dz = -1; dz <= 1; dz++) {
        let nx = cx + dx, ny = cy + dy, nz = cz + dz;
        if (dims) { nx = (nx + dims[0]) % dims[0]; ny = (ny + dims[1]) % dims[1]; nz = (nz + dims[2]) % dims[2]; }
        neighborKeys.add(key(nx, ny, nz));
      }
      for (const neighborKey of neighborKeys) {
        const neighbors = cells.get(neighborKey);
        if (!neighbors) continue;
        for (const i of members) for (const j of neighbors) {
          if (j <= i) continue;
          const a = atoms[i].pos, b = atoms[j].pos;
          const d = pbcDistance(a.x-b.x, a.y-b.y, a.z-b.z, this.sys.box);
          if (d.dx*d.dx + d.dy*d.dy + d.dz*d.dz <= maxR2) pairs.push([i, j]);
        }
      }
    }
    this.pairs = pairs;
    this.reference = atoms.map(a => a.pos.copy());
    this.valid = true;
  }
  getPairs() { if (this.needsRebuild()) this.rebuild(); return this.pairs; }
}

function nextPowerOfTwo(value) { let n=1; while (n < value) n <<= 1; return n; }
function erfc(x) {
  const z=Math.abs(x), t=1/(1+0.5*z);
  const ans=t*Math.exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))));
  return x >= 0 ? ans : 2-ans;
}

function fft1d(re, im, inverse) {
  const n=re.length;
  for (let i=1,j=0;i<n;i++) { let bit=n>>1; for (; j&bit; bit>>=1) j^=bit; j^=bit; if (i<j) { [re[i],re[j]]=[re[j],re[i]]; [im[i],im[j]]=[im[j],im[i]]; } }
  for (let len=2;len<=n;len<<=1) {
    const angle=(inverse ? 2 : -2)*Math.PI/len, wrStep=Math.cos(angle), wiStep=Math.sin(angle);
    for (let start=0;start<n;start+=len) { let wr=1, wi=0; const half=len>>1;
      for (let k=0;k<half;k++) { const even=start+k, odd=even+half, tr=re[odd]*wr-im[odd]*wi, ti=re[odd]*wi+im[odd]*wr; re[odd]=re[even]-tr; im[odd]=im[even]-ti; re[even]+=tr; im[even]+=ti; const nw=wr*wrStep-wi*wiStep; wi=wr*wiStep+wi*wrStep; wr=nw; }
    }
  }
  if (inverse) for (let i=0;i<n;i++) { re[i]/=n; im[i]/=n; }
}

class PMESolver {
  constructor(box) {
    this.box=box;
    this.nx=nextPowerOfTwo(Math.max(8, Math.ceil(box.x/config.pmeGridSpacing)));
    this.ny=nextPowerOfTwo(Math.max(8, Math.ceil(box.y/config.pmeGridSpacing)));
    this.nz=nextPowerOfTwo(Math.max(8, Math.ceil(box.z/config.pmeGridSpacing)));
    this.size=this.nx*this.ny*this.nz;
  }
  index(x,y,z) { return (x*this.ny+y)*this.nz+z; }
  transform(re, im, inverse) {
    const runLines=(count, lineLength, getter) => { const lr=new Float64Array(lineLength), li=new Float64Array(lineLength); for(let line=0;line<count;line++){ for(let k=0;k<lineLength;k++){const idx=getter(line,k);lr[k]=re[idx];li[k]=im[idx];} fft1d(lr,li,inverse); for(let k=0;k<lineLength;k++){const idx=getter(line,k);re[idx]=lr[k];im[idx]=li[k];} } };
    runLines(this.nx*this.ny,this.nz,(line,k)=>line*this.nz+k);
    runLines(this.nx*this.nz,this.ny,(line,k)=>this.index(Math.floor(line/this.nz),k,line%this.nz));
    runLines(this.ny*this.nz,this.nx,(line,k)=>this.index(k,Math.floor(line/this.nz),line%this.nz));
  }
  compute(atoms) {
    const rhoR=new Float64Array(this.size), rhoI=new Float64Array(this.size);
    const weights=[];
    for (const atom of atoms) {
      const gx=atom.pos.x/this.box.x*this.nx, gy=atom.pos.y/this.box.y*this.ny, gz=atom.pos.z/this.box.z*this.nz;
      const ix=Math.floor(gx), iy=Math.floor(gy), iz=Math.floor(gz), tx=gx-ix, ty=gy-iy, tz=gz-iz;
      const sites=[];
      for(let dx=0;dx<=1;dx++) for(let dy=0;dy<=1;dy++) for(let dz=0;dz<=1;dz++) { const w=(dx?tx:1-tx)*(dy?ty:1-ty)*(dz?tz:1-tz); const idx=this.index((ix+dx)%this.nx,(iy+dy)%this.ny,(iz+dz)%this.nz); rhoR[idx]+=atom.q*w; sites.push([idx,w]); }
      weights.push(sites);
    }
    this.transform(rhoR,rhoI,false);
    const phiR=new Float64Array(this.size), phiI=new Float64Array(this.size), exR=new Float64Array(this.size), exI=new Float64Array(this.size), eyR=new Float64Array(this.size), eyI=new Float64Array(this.size), ezR=new Float64Array(this.size), ezI=new Float64Array(this.size);
    const volume=this.box.x*this.box.y*this.box.z, alpha=config.ewaldAlpha, mesh=this.size;
    for(let x=0;x<this.nx;x++) for(let y=0;y<this.ny;y++) for(let z=0;z<this.nz;z++) { const idx=this.index(x,y,z); const kx=2*Math.PI*(x<=this.nx/2?x:x-this.nx)/this.box.x, ky=2*Math.PI*(y<=this.ny/2?y:y-this.ny)/this.box.y, kz=2*Math.PI*(z<=this.nz/2?z:z-this.nz)/this.box.z, k2=kx*kx+ky*ky+kz*kz; if(k2===0) continue; const g=mesh*4*Math.PI/volume*Math.exp(-k2/(4*alpha*alpha))/k2; phiR[idx]=g*rhoR[idx]; phiI[idx]=g*rhoI[idx]; exR[idx]=kx*phiI[idx]; exI[idx]=-kx*phiR[idx]; eyR[idx]=ky*phiI[idx]; eyI[idx]=-ky*phiR[idx]; ezR[idx]=kz*phiI[idx]; ezI[idx]=-kz*phiR[idx]; }
    this.transform(phiR,phiI,true); this.transform(exR,exI,true); this.transform(eyR,eyI,true); this.transform(ezR,ezI,true);
    const forces=atoms.map(()=>new Vec3()); let energy=0, q2=0;
    for(let i=0;i<atoms.length;i++) { let phi=0, ex=0, ey=0, ez=0; for(const [idx,w] of weights[i]) {phi+=w*phiR[idx];ex+=w*exR[idx];ey+=w*eyR[idx];ez+=w*ezR[idx];} const q=atoms[i].q; energy+=0.5*q*phi; q2+=q*q; forces[i].x=config.coulomb*q*ex; forces[i].y=config.coulomb*q*ey; forces[i].z=config.coulomb*q*ez; }
    return { forces, energy:config.coulomb*(energy-alpha*q2/Math.sqrt(Math.PI)) };
  }
}

function pbcDistance(dx, dy, dz, box) {
  if (!box) return { dx, dy, dz };
  dx -= Math.round(dx / box.x) * box.x;
  dy -= Math.round(dy / box.y) * box.y;
  dz -= Math.round(dz / box.z) * box.z;
  return { dx, dy, dz };
}

function applyPBC(sys) {
  if (!sys.box) return;
  for (const a of sys.atoms) {
    a.pos.x -= Math.floor(a.pos.x / sys.box.x) * sys.box.x;
    a.pos.y -= Math.floor(a.pos.y / sys.box.y) * sys.box.y;
    a.pos.z -= Math.floor(a.pos.z / sys.box.z) * sys.box.z;
  }
}

function randGaussian() {
  let u = 0, v = 0;
  while (u === 0) u = Math.random();
  while (v === 0) v = Math.random();
  return Math.sqrt(-2.0*Math.log(u))*Math.cos(2.0*Math.PI*v);
}

function constrainHydrogenBonds(sys) {
  const constrained = sys.bonds.filter(b => sys.atoms[b.i].element === 'H' || sys.atoms[b.j].element === 'H');
  for (let iteration = 0; iteration < config.constraintIterations; iteration++) {
    let largest = 0;
    for (const b of constrained) {
      const ai=sys.atoms[b.i], aj=sys.atoms[b.j];
      const d=pbcDistance(ai.pos.x-aj.pos.x, ai.pos.y-aj.pos.y, ai.pos.z-aj.pos.z, sys.box);
      const r2=d.dx*d.dx+d.dy*d.dy+d.dz*d.dz, target2=b.r0*b.r0, error=r2-target2;
      largest=Math.max(largest, Math.abs(error));
      if (Math.abs(error) < config.constraintTolerance) continue;
      const invMass=1/ai.mass+1/aj.mass, lambda=-error/(2*invMass*(r2+1e-16));
      ai.pos.x += lambda*d.dx/ai.mass; ai.pos.y += lambda*d.dy/ai.mass; ai.pos.z += lambda*d.dz/ai.mass;
      aj.pos.x -= lambda*d.dx/aj.mass; aj.pos.y -= lambda*d.dy/aj.mass; aj.pos.z -= lambda*d.dz/aj.mass;
    }
    if (largest < config.constraintTolerance) return;
  }
  throw new Error('SHAKE failed to converge; reduce dt or minimize the starting structure.');
}

function rattleHydrogenBonds(sys) {
  for (const b of sys.bonds) {
    const ai=sys.atoms[b.i], aj=sys.atoms[b.j];
    if (ai.element !== 'H' && aj.element !== 'H') continue;
    const d=pbcDistance(ai.pos.x-aj.pos.x, ai.pos.y-aj.pos.y, ai.pos.z-aj.pos.z, sys.box);
    const dvx=ai.vel.x-aj.vel.x, dvy=ai.vel.y-aj.vel.y, dvz=ai.vel.z-aj.vel.z;
    const lambda=(d.dx*dvx+d.dy*dvy+d.dz*dvz)/((d.dx*d.dx+d.dy*d.dy+d.dz*d.dz)*(1/ai.mass+1/aj.mass));
    ai.vel.x-=lambda*d.dx/ai.mass; ai.vel.y-=lambda*d.dy/ai.mass; ai.vel.z-=lambda*d.dz/ai.mass;
    aj.vel.x+=lambda*d.dx/aj.mass; aj.vel.y+=lambda*d.dy/aj.mass; aj.vel.z+=lambda*d.dz/aj.mass;
  }
}

function dihedralEnergy(dih, atoms) {
  const a = atoms[dih.i].pos, b = atoms[dih.j].pos, c = atoms[dih.k].pos, d = atoms[dih.l].pos;
  const ab = pbcDistance(b.x-a.x, b.y-a.y, b.z-a.z, null);
  const bc = pbcDistance(c.x-b.x, c.y-b.y, c.z-b.z, null);
  const cd = pbcDistance(d.x-c.x, d.y-c.y, d.z-c.z, null);
  const n1x = ab.dy*bc.dz-ab.dz*bc.dy, n1y = ab.dz*bc.dx-ab.dx*bc.dz, n1z = ab.dx*bc.dy-ab.dy*bc.dx;
  const n2x = bc.dy*cd.dz-bc.dz*cd.dy, n2y = bc.dz*cd.dx-bc.dx*cd.dz, n2z = bc.dx*cd.dy-bc.dy*cd.dx;
  const n1 = Math.hypot(n1x,n1y,n1z), n2 = Math.hypot(n2x,n2y,n2z), bn = Math.hypot(bc.dx,bc.dy,bc.dz);
  if (n1 < 1e-10 || n2 < 1e-10 || bn < 1e-10) return 0;
  const x = (n1x*n2x+n1y*n2y+n1z*n2z)/(n1*n2);
  const mx = n1y*bc.dz-n1z*bc.dy, my = n1z*bc.dx-n1x*bc.dz, mz = n1x*bc.dy-n1y*bc.dx;
  const y = (mx*n2x+my*n2y+mz*n2z)/(n1*n2*bn);
  const phi = Math.atan2(y, Math.max(-1, Math.min(1, x)));
  return dih.k * (1 + Math.cos(dih.n * phi - dih.delta));
}

function initVelocities(atoms) {
  for (const a of atoms) {
    const std = Math.sqrt((config.kB*config.targetTemp*config.fta)/a.mass);
    a.vel.x = std*randGaussian();
    a.vel.y = std*randGaussian();
    a.vel.z = std*randGaussian();
  }
}

async function minimizeGeometry(sys, steps = 250) {
  console.log("Minimizing geometry with L-BFGS...");
  const n = sys.atoms.length * 3, memory = 8;
  const dot = (a, b) => a.reduce((s, v, i) => s + v*b[i], 0);
  const coordinates = () => sys.atoms.flatMap(a => [a.pos.x, a.pos.y, a.pos.z]);
  const setCoordinates = x => sys.atoms.forEach((a, i) => { a.pos.x=x[3*i]; a.pos.y=x[3*i+1]; a.pos.z=x[3*i+2]; });
  const gradient = f => f.flatMap(v => [-v.x, -v.y, -v.z]);
  let state = await cf(sys), x = coordinates(), g = gradient(state.forces), history = [];
  for (let iter = 0; iter < steps; iter++) {
    const gNorm = Math.sqrt(dot(g,g)/n);
    if (gNorm < 1e-3) break;
    let q = g.slice(), alpha = [];
    for (let h = history.length-1; h >= 0; h--) { alpha[h] = history[h].rho*dot(history[h].s,q); q=q.map((v,k)=>v-alpha[h]*history[h].y[k]); }
    const last = history[history.length-1];
    const scale = last ? dot(last.s,last.y)/dot(last.y,last.y) : 1;
    let direction = q.map(v => -scale*v);
    for (let h = 0; h < history.length; h++) { const beta=history[h].rho*dot(history[h].y,direction); direction=direction.map((v,k)=>v+history[h].s[k]*(alpha[h]-beta)); }
    if (dot(direction,g) >= 0) direction = g.map(v => -v);
    const baseEnergy = state.energy.total, slope = dot(g,direction);
    let step = 0.02, next = null, nextX = null;
    for (let ls = 0; ls < 20; ls++) {
      nextX = x.map((v,k) => v + step*direction[k]); setCoordinates(nextX); applyPBC(sys); nextX = coordinates();
      next = await cf(sys);
      if (next.energy.total <= baseEnergy + 1e-4*step*slope) break;
      step *= 0.5;
    }
    if (!next || next.energy.total > baseEnergy) { setCoordinates(x); break; }
    const nextG = gradient(next.forces), s = nextX.map((v,k)=>v-x[k]), y=nextG.map((v,k)=>v-g[k]), sy=dot(s,y);
    if (sy > 1e-10) { history.push({s,y,rho:1/sy}); if (history.length > memory) history.shift(); }
    x=nextX; g=nextG; state=next;
  }
}

async function cf(sys) {
  const forces = sys.forces;
  for (let i = 0; i < forces.length; i++) {
    forces[i].x = 0; forces[i].y = 0; forces[i].z = 0;
  }
  const energy = { bond: 0, angle: 0, dihedral: 0, lj: 0, coulomb: 0, total: 0 };

  for (const b of sys.bonds) {
    const a1 = sys.atoms[b.i];
    const a2 = sys.atoms[b.j];
    let delta = pbcDistance(a1.pos.x-a2.pos.x, a1.pos.y-a2.pos.y, a1.pos.z-a2.pos.z, sys.box);
    const dx = delta.dx;
    const dy = delta.dy;
    const dz = delta.dz;
    const r = Math.sqrt(dx*dx+dy*dy+dz*dz)+1e-9;
    const stretch = r-b.r0;
    const s = (-b.k*stretch)/r;
    forces[b.i].x += dx*s;
    forces[b.i].y += dy*s;
    forces[b.i].z += dz*s;
    forces[b.j].x -= dx*s;
    forces[b.j].y -= dy*s;
    forces[b.j].z -= dz*s;
    energy.bond += 0.5*b.k*stretch*stretch;
  }

  for (const ang of sys.angles) {
    const a1 = sys.atoms[ang.i];
    const a2 = sys.atoms[ang.j];
    const a3 = sys.atoms[ang.k];
    let d1 = pbcDistance(a1.pos.x-a2.pos.x, a1.pos.y-a2.pos.y, a1.pos.z-a2.pos.z, sys.box);
    let d2 = pbcDistance(a3.pos.x-a2.pos.x, a3.pos.y-a2.pos.y, a3.pos.z-a2.pos.z, sys.box);
    const dx1 = d1.dx; const dy1 = d1.dy; const dz1 = d1.dz;
    const r1 = Math.sqrt(dx1*dx1+dy1*dy1+dz1*dz1)+1e-9;
    const dx2 = d2.dx; const dy2 = d2.dy; const dz2 = d2.dz;
    const r2 = Math.sqrt(dx2*dx2+dy2*dy2+dz2*dz2)+1e-9;
    let cosTheta = (dx1*dx2+dy1*dy2+dz1*dz2)/(r1*r2);
    cosTheta = Math.max(-1.0, Math.min(1.0, cosTheta));
    const theta = Math.acos(cosTheta);
    const dTheta = theta-ang.theta0;
    energy.angle += 0.5*ang.kTheta*dTheta*dTheta;
    let sinTheta = Math.sin(theta);
    if (Math.abs(sinTheta) < 1e-6) sinTheta = 1e-6;
    const factor = (-ang.kTheta*dTheta)/sinTheta;
    const f1x = factor*(cosTheta*(dx1/r1)-(dx2/r2))/r1;
    const f1y = factor*(cosTheta*(dy1/r1)-(dy2/r2))/r1;
    const f1z = factor*(cosTheta*(dz1/r1)-(dz2/r2))/r1;
    const f3x = factor*(cosTheta*(dx2/r2)-(dx1/r1))/r2;
    const f3y = factor*(cosTheta*(dy2/r2)-(dy1/r1))/r2;
    const f3z = factor*(cosTheta*(dz2/r2)-(dz1/r1))/r2;
    forces[ang.i].x += f1x;
    forces[ang.i].y += f1y;
    forces[ang.i].z += f1z;
    forces[ang.k].x += f3x;
    forces[ang.k].y += f3y;
    forces[ang.k].z += f3z;
    forces[ang.j].x -= (f1x+f3x);
    forces[ang.j].y -= (f1y+f3y);
    forces[ang.j].z -= (f1z+f3z);
  }

  for (const dih of sys.dihedrals) {
    energy.dihedral += dihedralEnergy(dih, sys.atoms);
    const h = 1e-5;
    for (const atomIndex of [dih.i, dih.j, dih.k, dih.l]) {
      const p = sys.atoms[atomIndex].pos;
      for (const axis of ['x', 'y', 'z']) {
        p[axis] += h; const plus = dihedralEnergy(dih, sys.atoms);
        p[axis] -= 2*h; const minus = dihedralEnergy(dih, sys.atoms);
        p[axis] += h;
        forces[atomIndex][axis] -= (plus - minus) / (2*h);
      }
    }
  }

  for (const [i, j] of sys.neighborList.getPairs()) {
    const pairClass = sys.pairClasses.get(`${i}:${j}`) || 0;
    const ai = sys.atoms[i], aj = sys.atoms[j];
    const d = pbcDistance(ai.pos.x-aj.pos.x, ai.pos.y-aj.pos.y, ai.pos.z-aj.pos.z, sys.box);
    const r2 = d.dx*d.dx+d.dy*d.dy+d.dz*d.dz;
    if (r2 >= config.cutoff*config.cutoff || r2 < 1e-12) continue;
    const r = Math.sqrt(r2);
    const t = Math.max(0, Math.min(1, (r-config.switchDistance)/(config.cutoff-config.switchDistance)));
    const sw = r <= config.switchDistance ? 1 : 1 - (10*t**3 - 15*t**4 + 6*t**5);
    const dSw = r <= config.switchDistance ? 0 : -(30*t*t - 60*t**3 + 30*t**4)/(config.cutoff-config.switchDistance);
    const li = ff.lj[ai.type] || ff.lj[ai.element] || ff.lj.undefDef;
    const lj = ff.lj[aj.type] || ff.lj[aj.element] || ff.lj.undefDef;
    const sig = 0.5*(li.sig+lj.sig), eps = Math.sqrt(li.eps*lj.eps);
    const sr2 = (sig*sig)/r2, sr6 = sr2**3, sr12 = sr6*sr6;
    const ljScale = pairClass === 2 ? config.oneFourLJScale : 1;
    const coulScale = pairClass === 2 ? config.oneFourCoulombScale : 1;
    const rawLJ = pairClass === 1 ? 0 : 4*eps*(sr12-sr6)*ljScale;
    const rawLJForce = pairClass === 1 ? 0 : (24*eps*(2*sr12-sr6)/r2)*ljScale;
    const qq=config.coulomb*ai.q*aj.q;
    const erfcTerm=erfc(config.ewaldAlpha*r), ewaldReal=qq*erfcTerm/r;
    const ewaldRealForce=qq*(erfcTerm/(r2*r)+(2*config.ewaldAlpha/Math.sqrt(Math.PI))*Math.exp(-((config.ewaldAlpha*r)**2))/r2);
    const rawCoul = sys.pme ? (pairClass === 0 ? ewaldReal : 0) : (pairClass === 1 ? 0 : qq/r*coulScale);
    const rawCoulForce = sys.pme ? (pairClass === 0 ? ewaldRealForce : 0) : (pairClass === 1 ? 0 : qq/(r2*r)*coulScale);
    const forceOverR = sw*rawLJForce - rawLJ*dSw/r + (sys.pme ? rawCoulForce : sw*rawCoulForce-rawCoul*dSw/r);
    const fx = d.dx*forceOverR, fy = d.dy*forceOverR, fz = d.dz*forceOverR;
    forces[i].x += fx; forces[i].y += fy; forces[i].z += fz;
    forces[j].x -= fx; forces[j].y -= fy; forces[j].z -= fz;
    energy.lj += sw*rawLJ;
    energy.coulomb += sys.pme ? rawCoul : sw*rawCoul;
  }

  if (sys.pme) {
    const reciprocal=sys.pme.compute(sys.atoms);
    energy.coulomb += reciprocal.energy;
    for (let i=0;i<sys.atoms.length;i++) forces[i].add(reciprocal.forces[i]);
    for (const [key, pairClass] of sys.pairClasses) {
      const [i,j]=key.split(':').map(Number), ai=sys.atoms[i], aj=sys.atoms[j];
      const d=pbcDistance(ai.pos.x-aj.pos.x,ai.pos.y-aj.pos.y,ai.pos.z-aj.pos.z,sys.box), r2=d.dx*d.dx+d.dy*d.dy+d.dz*d.dz;
      if(r2 < 1e-12) continue;
      const r=Math.sqrt(r2), qq=config.coulomb*ai.q*aj.q, erfTerm=1-erfc(config.ewaldAlpha*r);
      const fErf=qq*(erfTerm/(r2*r)-(2*config.ewaldAlpha/Math.sqrt(Math.PI))*Math.exp(-((config.ewaldAlpha*r)**2))/r2);
      const scale=pairClass===1 ? 0 : config.oneFourCoulombScale;
      const correction=qq*(scale/r-erfTerm/r), fCorrection=scale*qq/(r2*r)-fErf;
      const fx=d.dx*fCorrection, fy=d.dy*fCorrection, fz=d.dz*fCorrection;
      forces[i].x+=fx;forces[i].y+=fy;forces[i].z+=fz;forces[j].x-=fx;forces[j].y-=fy;forces[j].z-=fz;
      energy.coulomb+=correction;
    }
  }

  energy.total = energy.bond+energy.angle+energy.dihedral+energy.lj+energy.coulomb;
  return { forces, energy };
}

async function stepBAOAB(sys, forces) {
  const dt = config.dt;
  const halfDt = 0.5*dt;
  const decay = Math.exp(-config.gamma*dt);
  for (let i = 0; i < sys.atoms.length; i++) {
    const a = sys.atoms[i];
    const mult = config.fta/a.mass*halfDt;
    a.vel.x += forces[i].x*mult;
    a.vel.y += forces[i].y*mult;
    a.vel.z += forces[i].z*mult;
    a.pos.x += a.vel.x*halfDt;
    a.pos.y += a.vel.y*halfDt;
    a.pos.z += a.vel.z*halfDt;
  }
  constrainHydrogenBonds(sys);
  for (const a of sys.atoms) {
    const std = Math.sqrt((1-decay*decay)*config.kB*config.targetTemp*config.fta/a.mass);
    a.vel.x = decay*a.vel.x+std*randGaussian();
    a.vel.y = decay*a.vel.y+std*randGaussian();
    a.vel.z = decay*a.vel.z+std*randGaussian();
    a.pos.x += a.vel.x*halfDt;
    a.pos.y += a.vel.y*halfDt;
    a.pos.z += a.vel.z*halfDt;
  }

  applyPBC(sys);
  constrainHydrogenBonds(sys);

  const next = await cf(sys);
  for (let i = 0; i < sys.atoms.length; i++) {
    const a = sys.atoms[i];
    const mult = config.fta/a.mass*halfDt;
    a.vel.x += next.forces[i].x*mult;
    a.vel.y += next.forces[i].y*mult;
    a.vel.z += next.forces[i].z*mult;
  }
  rattleHydrogenBonds(sys);
  return next;
}

function calcTemp(atoms) {
  let ke = 0;
  for (const a of atoms) {
    ke += 0.5*a.mass*a.vel.magSq();
  }
  return (2*(ke/config.fta))/(3*atoms.length*config.kB);
}

const elementSymbols = new Set(Object.keys(masses).filter(k => k !== 'undefDef'));
const proteinAtomTypes = {
  N:'N_AMIDE', H:'H_AMIDE', CA:'C_ALPHA', C:'C_CARBONYL', O:'O_CARBONYL', OXT:'O_CARBONYL',
  CB:'C_ALIPHATIC', CG:'C_ALIPHATIC', CG1:'C_ALIPHATIC', CG2:'C_ALIPHATIC', CD:'C_ALIPHATIC', CD1:'C_ALIPHATIC', CD2:'C_ALIPHATIC', CE:'C_ALIPHATIC', CE1:'C_ALIPHATIC', CE2:'C_ALIPHATIC', CE3:'C_ALIPHATIC', CZ:'C_AROMATIC',
  ND1:'N_AROMATIC', ND2:'N_AMIDE', NE:'N_AMINE', NE1:'N_AROMATIC', NE2:'N_AMIDE', NZ:'N_AMINE', NH1:'N_AMINE', NH2:'N_AMINE',
  OG:'O_ALCOHOL', OG1:'O_ALCOHOL', OD1:'O_CARBONYL', OD2:'O_CARBONYL', OE1:'O_CARBONYL', OE2:'O_CARBONYL', OH:'O_ALCOHOL', SG:'S_THIOL', SD:'S_THIOETHER'
};

Object.assign(ff.lj, {
  N_AMIDE: {eps:0.17, sig:3.25}, N_AMINE:{eps:0.17, sig:3.25}, N_AROMATIC:{eps:0.17, sig:3.25},
  H_AMIDE: {eps:0.015, sig:1.10}, C_ALPHA:{eps:0.066, sig:3.50}, C_ALIPHATIC:{eps:0.066, sig:3.50}, C_AROMATIC:{eps:0.070, sig:3.55}, C_CARBONYL:{eps:0.066, sig:3.50},
  O_CARBONYL:{eps:0.1521, sig:3.15061}, O_ALCOHOL:{eps:0.1521, sig:3.15061}, S_THIOL:{eps:0.25, sig:3.60}, S_THIOETHER:{eps:0.25, sig:3.60}
});
const atomCharges = { N:-0.47, H:0.31, CA:0.07, HA:0.09, C:0.51, O:-0.51, OXT:-0.51, NZ:-0.30, NH1:-0.30, NH2:-0.30, NE:-0.30, ND2:-0.50, NE2:-0.50, OD1:-0.55, OD2:-0.55, OE1:-0.55, OE2:-0.55, OG:-0.42, OG1:-0.42, OH:-0.42, SG:-0.20, SD:-0.10 };
function inferElement(atomName, rawElement, record) {
  const specified = rawElement.trim();
  if (specified) { const e=specified[0].toUpperCase()+specified.slice(1).toLowerCase(); if (elementSymbols.has(e)) return e; }
  const letters = atomName.replace(/[^A-Za-z]/g, '');
  // PDB alignment distinguishes protein " CA " (carbon alpha) from "CA  " (calcium).
  if (record === 'ATOM' && proteinAtomTypes[atomName]) return atomName[0] === 'H' ? 'H' : ({N:'N',O:'O',S:'S'}[atomName[0]] || 'C');
  const two = letters.slice(0,2); const candidate = two[0]?.toUpperCase()+two.slice(1).toLowerCase();
  if (elementSymbols.has(candidate)) return candidate;
  const one = letters[0]?.toUpperCase();
  return elementSymbols.has(one) ? one : 'C';
}
function inferAtomType(atomName, resName, element) { return proteinAtomTypes[atomName] || `${resName}:${atomName}` || element; }
function inferCharge(atomName, resName, element) {
  if (atomCharges[atomName] !== undefined) return atomCharges[atomName];
  if (atomName.startsWith('H')) return 0.09;
  if (element === 'O') return -0.35;
  if (element === 'N') return -0.20;
  return 0.0;
}

function loadPDB(filePath) {
  const content = fs.readFileSync(filePath, 'utf-8');
  const lines = content.split('\n');
  const atoms = [];
  const bonds = [];
  const bondSet = new Set();
  const addBond = (i, j) => {
    const min = Math.min(i, j);
    const max = Math.max(i, j);
    const key = `${min}-${max}`;
    if (min !== max && !bondSet.has(key)) {
      bondSet.add(key);
      bonds.push([min, max]);
    }
  };
  const serialToIndex = new Map();
  const resGroups = [];
  let currentGroup = null;
  let box = null;

  for (const line of lines) {
    if (line.startsWith('CRYST1')) {
      const a = parseFloat(line.substring(6, 15).trim());
      const b = parseFloat(line.substring(15, 24).trim());
      const c = parseFloat(line.substring(24, 33).trim());
      if (!isNaN(a) && !isNaN(b) && !isNaN(c) && a > 0 && b > 0 && c > 0) {
        box = new Vec3(a, b, c);
      }
    } else if (line.startsWith('ATOM') || line.startsWith('HETATM')) {
      const record = line.substring(0, 6).trim();
      const serial = parseInt(line.substring(6, 11).trim(), 10);
      const atomName = line.substring(12, 16).trim();
      const resName = line.substring(17, 20).trim();
      const chainID = line.substring(21, 22).trim();
      const resSeq = parseInt(line.substring(22, 26).trim(), 10);
      const element = inferElement(atomName, line.substring(76, 78), record);
      const x = parseFloat(line.substring(30, 38).trim());
      const y = parseFloat(line.substring(38, 46).trim());
      const z = parseFloat(line.substring(46, 54).trim());
      const q = inferCharge(atomName, resName, element);
      const atom = new Atom(element, new Vec3(x, y, z), q, serial, inferAtomType(atomName, resName, element));
      const idx = atoms.length;
      serialToIndex.set(serial, idx);
      atoms.push(atom);
      if (!currentGroup || currentGroup.chainID !== chainID || currentGroup.resSeq !== resSeq) {
        currentGroup = { chainID, resSeq, resName, atomIdx: {} };
        resGroups.push(currentGroup);
      }
      currentGroup.atomIdx[atomName] = idx;
    } else if (line.startsWith('CONECT')) {
      const parts = line.substring(6).trim().split(/\s+/).map(v => parseInt(v, 10)).filter(v => !isNaN(v));
      if (parts.length > 1) {
        const srcSerial = parts[0];
        const srcIdx = serialToIndex.get(srcSerial);
        if (srcIdx !== undefined) {
          for (let i = 1; i < parts.length; i++) {
            const tgtIdx = serialToIndex.get(parts[i]);
            if (tgtIdx !== undefined) {
              addBond(srcIdx, tgtIdx);
            }
          }
        }
      }
    }
  }

  for (const group of resGroups) {
    for (const [a1, a2] of backboneBondPairs) {
      const i = group.atomIdx[a1];
      const j = group.atomIdx[a2];
      if (i !== undefined && j !== undefined) addBond(i, j);
    }
    const template = residueBondPairs[group.resName];
    if (template) {
      for (const [a1, a2] of template) {
        const i = group.atomIdx[a1];
        const j = group.atomIdx[a2];
        if (i !== undefined && j !== undefined) addBond(i, j);
      }
    }
    const names = Object.keys(group.atomIdx);
    const heavyNames = names.filter(n => atoms[group.atomIdx[n]].element !== 'H');
    const hNames = names.filter(n => atoms[group.atomIdx[n]].element === 'H');
    for (const hName of hNames) {
      const hBase = hName.slice(1);
      let bestName = null;
      let bestLen = -1;
      if (heavyNames.length === 1) {
        bestName = heavyNames[0];
      } else {
        for (const heavyName of heavyNames) {
          const heavyBase = heavyName.slice(1);
          if (heavyBase.length > 0 && hBase.startsWith(heavyBase) && heavyBase.length > bestLen) {
            bestName = heavyName;
            bestLen = heavyBase.length;
          }
        }
        if (bestName === null && heavyNames.includes('N')) {
          bestName = 'N';
        }
      }
      if (bestName !== null) {
        addBond(group.atomIdx[hName], group.atomIdx[bestName]);
      } else {
        const h = atoms[group.atomIdx[hName]];
        let nearest = null, nearestR2 = 1.25 * 1.25;
        for (const heavyName of heavyNames) {
          const heavy = atoms[group.atomIdx[heavyName]];
          const dx=h.pos.x-heavy.pos.x, dy=h.pos.y-heavy.pos.y, dz=h.pos.z-heavy.pos.z;
          const r2=dx*dx+dy*dy+dz*dz;
          if (r2 < nearestR2) { nearest=heavyName; nearestR2=r2; }
        }
        if (nearest !== null) addBond(group.atomIdx[hName], group.atomIdx[nearest]);
      }
    }
  }
  for (let g = 1; g < resGroups.length; g++) {
    const prev = resGroups[g-1];
    const curr = resGroups[g];
    if (prev.chainID === curr.chainID && curr.resSeq === prev.resSeq+1) {
      const i = prev.atomIdx['C'];
      const j = curr.atomIdx['N'];
      if (i !== undefined && j !== undefined) addBond(i, j);
    }
  }

  return { atoms, bonds, box };
}

function writeTrajectoryFrame(fd, atoms, step) {
  let frameData = `${atoms.length}\nFrame ${step}\n`;
  for (const a of atoms) {
    frameData += `${a.element.padEnd(2)} ${a.pos.x.toFixed(5)} ${a.pos.y.toFixed(5)} ${a.pos.z.toFixed(5)}\n`;
  }
  fs.writeSync(fd, frameData);
}

async function run(pdbPath, totalSteps, logInterval) {
  if (!fs.existsSync(pdbPath)) {
    console.error(`Error: PDB file not found at "${pdbPath}"`);
    return;
  }
  console.log(`\nLoading structure from ${pdbPath}...`);
  const { atoms, bonds, box } = loadPDB(pdbPath);
  console.log(`Loaded ${atoms.length} atoms.`);

  initVelocities(atoms);
  console.log(`Loaded ${bonds.length} bonds from topology/PDB.`);
  const sys = new System(atoms, bonds, box);
  await minimizeGeometry(sys, 250);
  const trajectoryFile = fs.openSync('traj.xyz', 'w');
  let { forces } = await cf(sys);
  const startTime = Date.now();

  function calcRg(atoms) {
    let totalMass = 0;
    let cx = 0, cy = 0, cz = 0;
    for (let i = 0; i < atoms.length; i++) {
      const a = atoms[i];
      totalMass += a.mass;
      cx += a.pos.x*a.mass;
      cy += a.pos.y*a.mass;
      cz += a.pos.z*a.mass;
    }
    cx /= totalMass; cy /= totalMass; cz /= totalMass;

    let sumSq = 0;
    for (let i = 0; i < atoms.length; i++) {
      const a = atoms[i];
      const dx = a.pos.x-cx;
      const dy = a.pos.y-cy;
      const dz = a.pos.z-cz;
      sumSq += a.mass*(dx*dx+dy*dy+dz*dz);
    }
    return Math.sqrt(sumSq/totalMass);
  }
  for (let s = 0; s < totalSteps; s++) {
    const res = await stepBAOAB(sys, forces);
    forces = res.forces;
    if (s % logInterval === 0) {
      writeTrajectoryFrame(trajectoryFile, sys.atoms, s);
      const elapsedSec = ((Date.now()-startTime)/1000).toFixed(1);
      const estRemainSec = s > 0 ? (((Date.now()-startTime)/s)*(totalSteps-s)/1000).toFixed(1) : "Calibrating...";
      const rg = calcRg(sys.atoms).toFixed(2);

      console.log(`Step ${s} | T: ${calcTemp(sys.atoms).toFixed(1)}K | E_tot: ${res.energy.total.toFixed(2)} | Rg: ${rg}Å | Elapsed: ${elapsedSec}s | ETA: ${estRemainSec}s`);
    }
  }

  fs.closeSync(trajectoryFile);
  console.log('\nSimulation complete. Trajectory saved to traj.xyz');
}

async function promptUser() {
  const rl = readline.createInterface({
    input: process.stdin,
    output: process.stdout
  });

  const ask = (query) => new Promise((resolve) => rl.question(query, resolve));

  const pdbPath = (await ask('PDB file name: [default: pdb/input.pdb]: ')).trim() || 'pdb/input.pdb';
  const stepsInput = (await ask('total simulation steps [default: 100]: ')).trim();
  const logInput = (await ask('logged frame per steps [default: 1]: ')).trim();

  rl.close();

  const totalSteps = parseInt(stepsInput, 10) || 100;
  const logInterval = parseInt(logInput, 10) || 1;

  run(pdbPath, totalSteps, logInterval);
}

if (require.main === module) promptUser();
module.exports = { Atom, Vec3, System, loadPDB, cf, minimizeGeometry, stepBAOAB, config };

/*

	Youve reached the end

*/
