(function (global, factory) {
	typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
	typeof define === 'function' && define.amd ? define(['exports'], factory) :
	(factory((global.j_lt = global.j_lt || {})));
}(this, (function (exports) { 'use strict';

var version="0.0.1"

function norm_pdf(mu, sigma) {
    var sqrt = Math.sqrt, pow = Math.pow, e = Math.E, pi = Math.PI;
    var a = 1 / (sigma * sqrt(2 * pi));
    return (function(xi) {
        return a*pow( e, - pow(xi - mu, 2) / (2 * pow(sigma, 2)) )
    });
}
function norm_max(sigma) { return(1/(sigma*Math.sqrt(2*Math.PI))) }

function gamma_pdf(mu,sigma) {
    var alpha=Math.pow(mu/sigma,2)
    var beta=mu/Math.pow(sigma,2)
    var k=Math.max(.0001,(alpha-1)/beta)
    var mx=Math.pow(k,alpha-1)*Math.pow(Math.E,-beta*k)
    return(function(xi){
        return Math.pow(xi,alpha-1)*Math.pow(Math.E,-beta*xi)/mx
    })
}

function lognormal_pdf(m,s) {
    const sigma=Math.sqrt(Math.log(s*s/(m*m)+1))
    const mu=Math.log(m)-sigma*sigma/2
    const k=Math.exp(mu-sigma*sigma)
    const mx=Math.exp(-Math.pow(Math.log(k)-mu,2)/(2*sigma*sigma))/(k*sigma)
    return(function(xi){
        return (1/mx)*Math.exp(-Math.pow(Math.log(xi)-mu,2)/(2*sigma*sigma))/(xi*sigma)
    })
}

function runSimSample(nSamp=50000){
    let r={}
    let noDelta={}
    let rMod={}
    for(var k in ioGroups){
        r[k]=d3.randomLogNormal(...aggregateParams_ln(ioGroups[k].baseline,false))
        noDelta[k]=function() {return(1)}
        rMod[k]=d3.randomLogNormal(...aggregateParams_ln(ioGroups[k].mods,false))
    }
    function est(mod) {
        let res=0
        for(k in ioGroups){
        res+=r[k]()*mod[k]()
        }
        return(res)
    }

    return([Array.from({length:nSamp},x=>est(noDelta)),Array.from({length:nSamp},x=>est(rMod))])
}

function runSim(nSamp=50000){ return(runSimApprox(nSamp)) }

function runSimApprox(){
      let r=[0,0]
      let rMod=[0,0]
      let a=[0,0],b=[0,0],ast=[0,0]
      for(var k in ioGroups){
        a=aggregateParams_ln(ioGroups[k].baseline,false)
        ast=getStats_ln(a[0],Math.pow(a[1],2))
        r[0]+=ast[0]
        r[1]+=Math.pow(ast[1],2)
        
        b=aggregateParams_ln(ioGroups[k].mods,false)
        b=getStats_ln(a[0]+b[0],Math.pow(a[1],2)+Math.pow(b[1],2))
        rMod[0]+=b[0]
        rMod[1]+=Math.pow(b[1],2)        
      }
    //   let smp=Array.from({length:nSamp},d3.randomNormal(r[0],Math.sqrt(r[1])))
    //   let delta=Array.from({length:nSamp},d3.randomNormal(rMod[0],Math.sqrt(rMod[1])))
      return([[r[0],Math.sqrt(r[1])],[rMod[0],Math.sqrt(rMod[1])]])
}

function updateElement(k,mk,gType){
    ioGroups[k][gType][mk].mu=Number(document.getElementById(gType+'-'+mk+'-mean-box').value)
    ioGroups[k][gType][mk].sigma=Number(document.getElementById(gType+'-'+mk+'-spread-box').value)
}

function kernelDensityEstimator(kernel, X) {
    return function(V) {
      return X.map(function(x) {
        return [x, d3.mean(V, function(v) { return kernel(x - v); }) ];
      });
    };
  }
function kernelEpanechnikov(k) {
    return function(v) {
        return Math.abs(v /= k) <= 1 ? 0.75 * (1 - v * v) / k : 0;
    };
}

function updateEstimates(sample=true) {
    // for(var k in ioGroups){
    //     for(var mk in ioGroups[k].mods){
    //         ioGroups[k].mods[mk].mu=Number(document.getElementById('delta-'+mk+'-box').value)
    //     }
    //     ioGroups[k].mods
    // }
    let den=null, modDen=null,data=null,deltaDen=null
    if(sample){
        data = runSimSample()
        outGrp.updateData(data)
        den=outGrp.density
        modDen=outGrp.modDensity
        ap.updateX(Math.min(Math.min(...den.map(x=>x[0])),Math.min(...modDen.map(x=>x[0])))
            ,Math.max(Math.max(...den.map(x=>x[0])),Math.max(...modDen.map(x=>x[0]))),0)
        ap.updateY(0,Math.max(Math.max(...den.map(x=>x[1])),Math.max(...modDen.map(x=>x[1]))),0)

        var delta=Array.from({length:data[0].length},x=>0)
        for(var i=0;i<delta.length;i++) delta[i]=data[1][i]-data[0][i]
        let kde=function(XX,k=10) { return(kernelDensityEstimator(kernelEpanechnikov(k), XX)) }
        let ff=kde(d3.scaleLinear().domain([Math.min(...delta),Math.max(...delta)]).ticks(100))
        deltaDen=ff(delta)
        let alpha=deltaDen.map(x=>x[1]).reduce((a,b)=>a+b)*(deltaDen[1][0]-deltaDen[0][0])
        deltaDen.map(x=>[x[0],x[1]/alpha])
        sim_delta.updateX(Math.min(...deltaDen.map(x=>x[0])),Math.max(...deltaDen.map(x=>x[0])),0)
        sim_delta.updateY(0,Math.max(...deltaDen.map(x=>x[1])),0)
    } else {
        let params=runSimApprox()
        den= d3.scaleLinear().domain([outGrp.xMin,outGrp.xMax])
            .ticks(100).map(x=>[x,norm_pdf(...params[0])(x)])
        modDen= d3.scaleLinear().domain([outGrp.xMin,outGrp.xMax])
            .ticks(100).map(x=>[x,norm_pdf(...params[1])(x)])
        let mu=params[1][0]-params[0][0]
        let sig=Math.sqrt(Math.pow(params[0][1],2)+Math.pow(params[1][1],2))
        deltaDen= d3.scaleLinear().domain([mu-6*sig,mu+6*sig]).ticks(100)
            .map(x=>[x,norm_pdf(mu,sig)(x)])
        ap.updateY(0,Math.max(norm_max(params[0][1]),norm_max(params[1][1])))
        sim_delta.updateY(0,norm_max(sig))
    }
    ap.curves[1].datum(den).attr("d",ap.area)
    ap.curves[0].datum(modDen).attr("d",ap.area)
    sim_delta.curves[0].datum(deltaDen).attr("d",sim_delta.area)
}
  
function selectCC(ccName){
    currentCC=ccName
    document.getElementById('blLabel').innerHTML=currentCC
    let loc=document.getElementById('baselineAssumptions')
    while(loc.firstChild) { loc.removeChild(loc.firstChild)}
    loc=document.getElementById('interventionImpacts')
    while(loc.firstChild) { loc.removeChild(loc.firstChild)}
    for(var mk in ioGroups[currentCC].baseline){
        j_lt.putInputGroup('baselineAssumptions','baseline-'+mk,ioGroups[currentCC].baseline[mk]
            ,function(){j_lt.updateElement(currentCC,mk,'baseline'); j_lt.updateEstimates()},0)
    }
    for(mk in ioGroups[currentCC].mods){
        j_lt.putInputGroup("interventionImpacts",'mods-'+mk,ioGroups[currentCC].mods[mk]
            ,function(){j_lt.updateElement(currentCC,mk,'mods'); j_lt.updateEstimates()},0)
            
    }
}

function putInputGroup(targetID,label,grp,upFunc,tDelay=500) {
    grp.updateDensity()
    grp.updateX()
    let ap=buildAreaPlot(targetID,label,1,100,460,grp.title)
    ap.updateX(grp.xMin,grp.xMax,tDelay)
    ap.curves[0].datum(grp.density).attr("d",ap.area)

    let upd=function(label,tDelay){
        grp.mu=Number(document.getElementById(label+'-mean-box').value)
        grp.sigma=Number(document.getElementById(label+'-spread-box').value)
        grp.updateDensity()
        ap.curves[0].datum(grp.density).transition().duration(tDelay).attr("d",ap.area)
    }
    let upx=function(label){
        grp.mu=Number(document.getElementById(label+'-mean-box').value)
        grp.sigma=Number(document.getElementById(label+'-spread-box').value)
        grp.updateX()
        ap.updateX(grp.xMin,grp.xMax)
        upFunc(label)
    }

    addGroupInput(targetID,label+'-mean','Mean',grp.mu,grp.xMin,grp.xMax
        ,function() {upx(label); upd(label,tDelay)},function(){upd(label,0)})
    addGroupInput(targetID,label+'-spread','Spread',grp.sigma,.0001,grp.sigMax
        ,function(){upx(label); upd(label,tDelay)},function(){upd(label,0)})
}

function aggregateParams_ln(grpList,returnStats=true) {
    let m=0
    let s=0
    for(var nm in grpList){
        let gg=grpList[nm]
        const params=getParams_ln(gg.mu,gg.sigma)
        // const si=Math.log(gg.sigma*gg.sigma/(gg.mu*gg.mu)+1)
        // m+=Math.log(gg.mu)-si/2
        // s+=si
        m+=params[0]
        s+=params[1]
    }
    // let mu=Math.exp(m+s/2)
    // let sigma=Math.sqrt((Math.exp(s)-1)*Math.exp(2*m+s))
    if(returnStats) {
        return(getStats_ln(m,s))
    } else {
        return([m,Math.sqrt(s)])
    }
    
}

function aggregateGroup(key,grpList) {
    let sts=aggregateParams_ln(grpList)
    let grp=grpConstructor(key.replaceAll(' ','_'),key,sts[0],sts[1],'lognormal')
    grp.updateX()
    return(grp)
}

function grpConstructor(label,title,mu,sigma=0,density='lognormal'){
    if(sigma==0) sigma=mu/100
    var grp={label:label,
        title:title,
        mu:mu, sigma:sigma,
        yMin:0,yMax:1,
        sigMax:sigma*6,density:density}
    if(density==="normal") {
        grp['pdf']=norm_pdf
        grp['rand']=function() {
            return d3.randomNormal(grp.mu,grp.sigma)()
        }
    } else if (density==='lognormal') {
        grp['pdf']=lognormal_pdf
        grp['rand']=function() {
            const s=Math.sqrt(Math.log(grp.sigma*grp.sigma/(grp.mu*grp.mu)+1))
            const m=Math.log(grp.mu)-s*s/2
            return d3.randomLogNormal(m,s)()
        }
    }
    grp.updateDensity=function() {
        let ff=grp.pdf(grp.mu,grp.sigma)
        grp.density= d3.scaleLinear().domain([grp.xMin,grp.xMax])
            .ticks(100).map(x=>[x,ff(x)])
    }
    grp.updateX=function() {
        grp.xMin=Math.max(0.0001,grp.mu-6*grp.sigma)
        grp.xMax=grp.mu+6*grp.sigma
    }
    grp.updateX()
    grp.updateDensity()
    return(grp)
}

function outputTag(label,title,data,plotFunc='densityPlot') {
    var grp={label:label,title:title}
    if(plotFunc==='densityPlot') {
        grp['plot']=function(XX,k=10) { return(kernelDensityEstimator(kernelEpanechnikov(k), XX)) }
    }
    grp.updateDensity= function() {
        let ff=grp.plot(d3.scaleLinear().domain([grp.xMin,grp.xMax]).ticks(100))
        grp.density=ff(grp.data)
        grp.modDensity=ff(grp.modData)
    }
    grp.updateData=function(data) {
        grp.data=data[0]
        grp.modData=data[1]
        // grp['xMax'] = data.sort(function(a, b){return b-a}).slice(nn,nn+1)[0]
        grp['xMax'] = Math.max(...grp.data.concat(grp.modData))
        grp['xMin']=Math.min(...grp.data.concat(grp.modData))
        grp.updateDensity()
        grp['yMax']=Math.max(Math.max(...grp.density.map(x=>x[1])),Math.max(...grp.modDensity.map(x=>x[1])))
    }
    grp['yMin']=0;
    grp.updateData(data)
    return(grp);
}

function getParams_ln(mn,sd){
    const s=Math.log(Math.pow(sd/mn,2)+1)
    return([Math.log(mn)-s/2,s])
}
function getStats_ln(alpha,phi){
    return([Math.exp(alpha+phi/2),Math.sqrt((Math.exp(phi)-1)*Math.exp(2*alpha+phi))])
}

function updateIOModGroup(kk,gType="mods",newMean=NaN) {
    let mp=[0,0]
    let n=0
    for(var v in ioGroups[kk][gType]){
        let nv=getParams_ln(ioGroups[kk][gType][v].mu,ioGroups[kk][gType][v].sigma)
        mp[0]+=nv[0]
        mp[1]+=nv[1]
        n++
    }
    let nw=getStats_ln(...mp)
    if(isNaN(newMean)) newMean=Number(document.getElementById(kk.replaceAll(' ','_')+'-box').value)
    let newSD=nw[1]*newMean/nw[0]
    let delta=Math.log(newMean/nw[0])
    let b=(Math.log(Math.exp(-2*delta)*Math.pow(newSD/nw[1],2)*(Math.exp(mp[1])-1)+1)-mp[1])/(2*n)
    let a=(delta-n*b)/n
    for(v in ioGroups[kk][gType]){
        nw=getParams_ln(ioGroups[kk][gType][v].mu,ioGroups[kk][gType][v].sigma)
        nw=getStats_ln(nw[0]+a,nw[1]+2*b)
        ioGroups[kk][gType][v].mu=nw[0]
        ioGroups[kk][gType][v].sigma=nw[1]
        ioGroups[kk][gType][v].updateX()
        ioGroups[kk][gType][v].updateDensity()
    }
}

function updateIOGroup(kk,gType='baseline') {
    let mp=[0,0]
    let n=0
    for(var v in ioGroups[kk][gType]){
        let nv=getParams_ln(ioGroups[kk][gType][v].mu,ioGroups[kk][gType][v].sigma)
        mp[0]+=nv[0]
        mp[1]+=nv[1]
        n++
    }
    let nw=getStats_ln(...mp)
    let newMean=Number(document.getElementById(kk.replaceAll(' ','_')+'-mean-box').value)
    let newSD=Number(document.getElementById(kk.replaceAll(' ','_')+'-spread-box').value)
    let delta=Math.log(newMean/nw[0])
    let b=(Math.log(Math.exp(-2*delta)*Math.pow(newSD/nw[1],2)*(Math.exp(mp[1])-1)+1)-mp[1])/(2*n)
    let a=(delta-n*b)/n
    for(v in ioGroups[kk][gType]){
        nw=getParams_ln(ioGroups[kk][gType][v].mu,ioGroups[kk][gType][v].sigma)
        nw=getStats_ln(nw[0]+a,nw[1]+2*b)
        ioGroups[kk][gType][v].mu=nw[0]
        ioGroups[kk][gType][v].sigma=nw[1]
        ioGroups[kk][gType][v].updateX()
        ioGroups[kk][gType][v].updateDensity()
    }
}
  
function addGroupInput(targetID,label,title,value,min,max,changeFunc=function(){},inputFunc=function(){}) {
    let mdv=document.createElement("div");
    mdv.setAttribute("class","inputGroup")
    let lbl=document.createElement("label")
    lbl.setAttribute("class","inputGroupLabel")
    lbl.appendChild(document.createTextNode(title))
    mdv.appendChild(lbl)
    let rng=document.createElement("input");
    rng.type='range'
    rng.setAttribute("class","inputGroupSlide")
    rng.id=label+'-slide';
    rng.min=min;
    rng.max=max
    if(rng.max-rng.min<100) { rng.step=(rng.max-rng.min)/1000; }
    rng.value=value;
    mdv.appendChild(rng);
    let bx=document.createElement('input');
    bx.setAttribute("class","inputGroupBox")
    bx.value=value;
    bx.id=label+'-box'
    mdv.appendChild(bx);
    document.getElementById(targetID).appendChild(mdv)
    d3.select("#"+label+"-box").on("change", function(){
        let sld=document.getElementById(label+"-slide")
        if(Number(this.value)>Number(sld.max)) { sld.max=this.value*1.1-sld.min*.1 }
        if(Number(this.value)<Number(sld.min)) { sld.min=this.value*1.1-sld.max*.1 }
        sld.value=Number(this.value)
        changeFunc(label)
    })
    d3.select("#"+label+"-slide").on("input", function(){
        document.getElementById(label+"-box").value=Number(this.value)
        inputFunc(label)
    })
    d3.select("#"+label+"-slide").on("change", function(){
        changeFunc(label)
    })
}

function buildAreaPlot(targetID,label,nCurves=1,hh=100,ww=460,title='') {
    var margin = {top: 30, right: 30, bottom: 30, left: 50},
    width = ww - margin.left - margin.right,
    height = hh - margin.top - margin.bottom;

    let svg=insertAreaPlot(targetID,label,margin,hh,ww,title)

    var x = d3.scaleLinear()
        .domain([0,1])
        .range([0, width]);
    svg.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + height + ")")
        .call(d3.axisBottom(x));

    var y = d3.scaleLinear()
        .domain([0,1])
        .range([height, 0]);
    svg.append("g")
        .attr("class", "y axis")
        .call(d3.axisLeft(y))
        .style("opacity",0)

    let colors=["#0f62a0","#800000","#B6B7B6","#30305f","#101010","#69b3a2"]
    let curves=[]
    for(var i=0;i<nCurves;i++) {
        curves.push(svg.append("g")
            .append("path")
            .attr('id','curve_'+i)
            .attr("fill", colors[i])
            .attr("opacity", ".5")
        //     .attr("stroke", "#000")
        //     .attr("stroke-width", 1)
        //     .attr("stroke-linejoin", "round")
        //     .attr("d",  d3.line()
        //     .curve(d3.curveBasis)
        //         .x(function(d) { return x(d[0]); })
        //         .y(function(d) { return y(d[1]); })
        )
    }

    var grp={label:label,height:height,width:width,svg:svg,x:x,y:y,curves:curves}
    grp.area = d3.area()
        .x(function(d) {return grp.x(d[0]); })
        .y0(function(d) {return grp.y(0)})
        .y1(function(d) { return grp.y(d[1]); })

    grp.updateY=function(mn,mx,tDelay=500) { 
        grp.y.domain([mn,mx])
        grp.svg.select(".y").transition().duration(tDelay)
            .call(d3.axisLeft(grp.y))
    }
    grp.updateX=function(mn,mx,tDelay=500) {
        grp.x.domain([mn,mx])
        grp.svg.select(".x").transition().duration(tDelay)
            .call(d3.axisBottom(grp.x))
    }
    return(grp)
}

function insertAreaPlot(targetID,label,margin,hh,ww,title) {
    var svgWindow=document.createElement('div');
    svgWindow.id=label+'-areaPlot';
    document.getElementById(targetID).appendChild(svgWindow)

// append the svg object to the body of the page
    var svg = d3.select("#"+label+"-areaPlot")
        .append("svg")
            .attr("width", ww)
            .attr("height", hh)
        .append("g")
            .attr("transform","translate(" + margin.left + "," + margin.top + ")");
        svg.append("text")
            .attr("x", ( (ww - margin.left - margin.right) / 2))             
            .attr("y", 0 - (margin.top / 2))
            .attr("text-anchor", "middle")  
            .style("font-size", "16px") 
            .text(title);
        return(svg)
}

function selectTab(tabName) {
    var tabcontent = document.getElementsByClassName("tabcontent");
    for (var i = 0; i < tabcontent.length; i++) {
      tabcontent[i].style.display = "none";
    }
    document.getElementById(tabName).style.display = "block";
    if(tabName==='Cost Categories') selectCC(currentCC)
    if(tabName==='Populations') {
        let loc=document.getElementById('patientDetailsAll')
        while(loc.firstChild) { loc.removeChild(loc.firstChild)}
        for(k in ioGroups){
            j_lt.putInputGroup('patientDetailsAll',k.replaceAll(' ','_')
            ,j_lt.aggregateGroup(k,ioGroups[k].baseline)
            ,function(k){j_lt.updateIOGroup(k.replaceAll('_',' '));j_lt.updateEstimates()},0)
        }
    }
    if(tabName==='Estimates'){
        for(k in ioGroups){
            let gg=j_lt.aggregateGroup(k,ioGroups[k].mods)
            document.getElementById(k.replaceAll(' ','_')+'-box').value=gg.mu
            document.getElementById(k.replaceAll(' ','_')+'-slide').value=gg.mu
            document.getElementById(k.replaceAll(' ','_')+'-slide').max=Math.max(gg.xMax,2)
        }
        // j_lt.addGroupInput("Estimates",k.replaceAll(' ','_'),k,1,0,2
        // ,function(k){j_lt.updateIOModGroup(k.replaceAll('_',' ')); j_lt.updateEstimates()})

    }
}

function selectIVN(){
    v=document.getElementById('ivnList').value
    for(let k in ioGroups){
        for(let mt in ioGroups[k].mods){
            ioGroups[k].mods[mt].mu=1
            ioGroups[k].mods[mt].sigma=.01
        }
        document.getElementById(k.replaceAll(' ','_')+'-box').value=1
        document.getElementById(k.replaceAll(' ','_')+'-slide').value=1
    }
    for(k in ivnList[v].modifiers){
      let mu=ivnList[v].modifiers[k]
      updateIOModGroup(k,"mods",mu)
      document.getElementById(k.replaceAll(' ','_')+'-box').value=mu
      document.getElementById(k.replaceAll(' ','_')+'-slide').value=mu
    }
    j_lt.updateEstimates()
  }  
  
exports.norm_pdf = norm_pdf;
exports.version=version;
exports.grpConstructor=grpConstructor;
exports.putInputGroup=putInputGroup;
exports.gamma_pdf=gamma_pdf;
exports.lognormal_pdf=lognormal_pdf;
exports.kernelDensityEstimator=kernelDensityEstimator;
exports.kernelEpanechnikov=kernelEpanechnikov;
exports.outputTag=outputTag;
exports.buildAreaPlot=buildAreaPlot;
exports.addGroupInput=addGroupInput;
exports.aggregateGroup=aggregateGroup;
exports.updateIOGroup=updateIOGroup;
exports.updateIOModGroup=updateIOModGroup;
exports.updateEstimates=updateEstimates;
exports.selectTab=selectTab;
exports.selectCC=selectCC;
exports.updateElement=updateElement;
exports.runSim=runSim;
exports.runSimSample=runSimSample;
exports.runSimApprox=runSimApprox;
exports.selectIVN=selectIVN;
exports.aggregateParams_ln=aggregateParams_ln;
exports.getStats_ln=getStats_ln;

Object.defineProperty(exports, '__esModule', { value: true });

})));