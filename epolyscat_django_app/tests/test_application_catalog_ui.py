import subprocess
from pathlib import Path


def test_frontend_intersects_discovery_with_supported_contracts():
    script = r'''
const fs = require("fs"), babel = require("@babel/core"), assert = require("assert");
const source = fs.readFileSync("src/utils/application-catalog.js", "utf8");
const compiled = babel.transformSync(source, {plugins:["@babel/plugin-transform-modules-commonjs"]}).code;
const m={exports:{}};new Function("module","exports",compiled)(m,m.exports);
const contracts=[
 {id:"module",applications:[{id:"ePolyScat",requiredFiles:["input"]},{id:"Gaussian16"}]},
 {id:"utility",applications:[{id:"CnvMath"},{id:"CnvLinFull"}]},
 {id:"workflow",applications:[{id:"Data_Gen",workflowApplicationIds:["Gaussian16"]},{id:"Analysis",workflowApplicationIds:["CnvMath","CnvLinFull"]}]}
];
const apps=[{id:"ePolyScat",kind:"module",label:"Registered ePolyScat",moduleId:"fresh-id"},{id:"CnvLinFull",kind:"utility"},{id:"Unknown",kind:"module"}];
const result=m.exports.filterApplicationContracts(contracts,apps);
assert.deepEqual(result[0].applications.map(a=>a.id),["ePolyScat"]);
assert.equal(result[0].applications[0].moduleId,"fresh-id");
assert.deepEqual(result[0].applications[0].requiredFiles,["input"]);
assert.deepEqual(result[1].applications.map(a=>a.id),["CnvLinFull"]);
assert.deepEqual(result[2].applications[0].workflowApplicationIds,[]);
assert.deepEqual(result[2].applications[1].workflowApplicationIds,["CnvLinFull"]);
assert.equal(m.exports.filterApplicationContracts(contracts,[])[0].applications.length,0);
assert.equal(contracts[0].applications.length,2);
assert.deepEqual(m.exports.applicationSelectionForRun({runMode:"module",moduleApplication:"RemovedGaussian"}),{runType:"module",applicationId:"RemovedGaussian"});
assert.deepEqual(m.exports.applicationSelectionForRun({runMode:"utility",utilityApplication:"MoldenMerge"}),{runType:"utility",applicationId:"MoldenMerge"});
assert.equal(m.exports.applicationSelectionForRun({runMode:"workflow",workflowStage:"Analysis",utilityApplication:"CnvLinFull"}).stageSelections.Analysis,"CnvLinFull");
'''
    result = subprocess.run(["node", "-e", script], cwd=Path(__file__).resolve().parents[1], text=True, capture_output=True)
    assert result.returncode == 0, result.stderr


def test_catalog_store_deduplicates_requests_and_clears_stale_choices_on_failure():
    script = r'''
const fs=require("fs"), assert=require("assert");
const source=fs.readFileSync("src/store/modules/settings.store.js","utf8").replace(/^import .*;$/gm, "").replace("export default", "return");
let resolve, calls=0;
const service={all:()=>{calls++;return new Promise(r=>{resolve=r;});}};
const store=new Function("SettingsService",source)(service);
const context={commit:(mutation,value)=>store.mutations[mutation](store.state,value)};
(async()=>{
 const first=store.actions.fetchSettings(context), second=store.actions.fetchSettings(context);
 assert.equal(calls,1); assert.equal(store.state.loading,true);
 resolve({applications:[{id:"ePolyScat",moduleId:"new-id"}]});
 assert.deepEqual(await Promise.all([first,second]),[true,true]);
 assert.equal(store.state.loaded,true);
 service.all=async()=>{throw Error("offline");};
 assert.equal(await store.actions.fetchSettings(context),false);
 assert.equal(store.state.loaded,false); assert.deepEqual(store.state.settings,{});
 assert.ok(store.state.error);
 service.all=async()=>({applications:[]});
 assert.equal(await store.actions.fetchSettings(context),true);
 assert.equal(store.state.error,""); assert.equal(store.state.loaded,true);
})().catch(error=>{console.error(error);process.exitCode=1;});
'''
    result = subprocess.run(["node", "-e", script], cwd=Path(__file__).resolve().parents[1], text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
