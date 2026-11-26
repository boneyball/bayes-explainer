// Verifies the page's analytic math against a Monte Carlo simulation.
const fs=require('fs');
const html=fs.readFileSync(__dirname+'/../index.html','utf8');
const code=html.split('/*MATH-START*/')[1].split('/*MATH-END*/')[0];
const M=new Function(code+';return {betaCdf,betaPdf,betaQuantile,probBBeatsA,liftSummary};')();

function rng(seed){return()=>{seed|=0;seed=seed+0x6D2B79F5|0;let t=Math.imul(seed^seed>>>15,1|seed);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;};}
const rand=rng(12345);
function normal(){let u=0;while(!u)u=rand();return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*rand());}
function gamma(k){ // Marsaglia-Tsang, k>=1
  const d=k-1/3,c=1/Math.sqrt(9*d);
  for(;;){let x,v;do{x=normal();v=1+c*x;}while(v<=0);v=v*v*v;const u=rand();
    if(Math.log(u)<0.5*x*x+d-d*v+d*Math.log(v))return d*v;}
}
const betaDraw=(a,b)=>{const x=gamma(a),y=gamma(b);return x/(x+y);};

const cases=[[1000,100,1000,118],[1000,100,1000,100],[200,20,200,31],[50,4,60,9],[20000,2000,20000,2100],[10,1,10,3]];
const N=2_000_000; let worst=0;
console.log('case (A v/c, B v/c) | P(B>A) analytic vs MC | lift mean | lift 95% CI analytic vs MC');
for(const [av,ac,bv,bc] of cases){
  const aA=1+ac,bA=1+av-ac,aB=1+bc,bB=1+bv-bc;
  const p=M.probBBeatsA(aA,bA,aB,bB), L=M.liftSummary(aA,bA,aB,bB);
  let wins=0,sum=0; const lifts=new Float64Array(N);
  for(let i=0;i<N;i++){const a=betaDraw(aA,bA),b=betaDraw(aB,bB);if(b>a)wins++;lifts[i]=b/a-1;sum+=lifts[i];}
  lifts.sort(); const mc={p:wins/N,mean:sum/N,lo:lifts[Math.floor(.025*N)],hi:lifts[Math.floor(.975*N)]};
  const se=Math.sqrt(mc.p*(1-mc.p)/N);
  worst=Math.max(worst,Math.abs(p-mc.p));
  console.log(`A ${ac}/${av} B ${bc}/${bv} | ${(p*100).toFixed(2)}% vs ${(mc.p*100).toFixed(2)}% (MC s.e. ${(se*100).toFixed(2)}) | mean ${(L.mean*100).toFixed(2)}% vs ${(mc.mean*100).toFixed(2)}% | CI [${(L.lo*100).toFixed(1)}, ${(L.hi*100).toFixed(1)}] vs [${(mc.lo*100).toFixed(1)}, ${(mc.hi*100).toFixed(1)}]`);
}
// CDF/quantile sanity: Beta(1,1) is uniform; Beta(2,5) median known ≈ 0.26445
console.log('uniform cdf(.3)=',M.betaCdf(.3,1,1),' Beta(2,5) median=',M.betaQuantile(.5,2,5).toFixed(5),'(expect 0.26445)');
console.log('worst |P(B>A) analytic - MC| =',(worst*100).toFixed(3),'points');
process.exit(worst<0.003?0:1);
