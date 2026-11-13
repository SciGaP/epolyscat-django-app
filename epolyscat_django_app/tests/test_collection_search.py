import subprocess
from pathlib import Path

from django.contrib.auth import get_user_model
from django.test import TestCase
from rest_framework.test import APIRequestFactory, force_authenticate

from epolyscat_django_app import models, views


class CollectionSearchTests(TestCase):
    def setUp(self):
        self.user = get_user_model().objects.create_user(username="search-owner")
        self.other = get_user_model().objects.create_user(username="search-other")

    def listing(self, viewset, **params):
        request = APIRequestFactory().get("/api/collection/", params)
        force_authenticate(request, user=self.user)
        result = viewset.as_view({"get": "list"})(request)
        self.assertEqual(result.status_code, 200)
        return result.data

    def test_experiment_search_precedes_pagination_and_respects_ownership(self):
        match = models.Experiment.objects.create(name="Target", description="Helium spectrum", owner=self.user)
        for i in range(12):
            models.Experiment.objects.create(name=f"Other {i}", description="unrelated", owner=self.user)
        models.Experiment.objects.create(name="Helium private", description="", owner=self.other)
        models.Experiment.objects.create(name="Helium deleted", description="", owner=self.user, deleted=True)
        for term in ("HELIUM", "Target"):
            result = self.listing(views.ExperimentViewSet, search=term, page_size=2)
            self.assertEqual(result["count"], 1)
            self.assertEqual([row["id"] for row in result["results"]], [match.pk])
        self.assertEqual(self.listing(views.ExperimentViewSet, search="missing")["count"], 0)
        self.assertEqual(self.listing(views.ExperimentViewSet, search="")["count"], 13)

    def test_view_search_filters_all_pages_and_preserves_owner_scope(self):
        for i in range(12):
            models.View.objects.create(name=f"Other {i}", owner=self.user)
        target = models.View.objects.create(name="Helium collection", owner=self.user)
        models.View.objects.create(name="Helium private", owner=self.other)
        result = self.listing(views.ViewsViewSet, search="HELIUM", page_size=2)
        self.assertEqual(result["count"], 1)
        self.assertEqual([row["id"] for row in result["results"]], [target.pk])
        self.assertEqual(self.listing(views.ViewsViewSet, search="missing")["count"], 0)
        self.assertEqual(self.listing(views.ViewsViewSet, search="")["count"], 13)


def test_search_results_and_totals_are_cached_per_query():
    script = r'''
const fs=require("fs"),assert=require("assert");
function load(file, serviceName, service) {
 const code=fs.readFileSync(file,"utf8").replace(/^import .*;$/gm,"").replace("export default","return");
 return new Function(serviceName,code)(service);
}
(async()=>{
 for(const kind of ["experiment","view"]) {
  const plural=kind==="experiment"?"Experiments":"Views", id=kind+"Id", requests=[];
  const service={ ["fetchAll"+plural]: args=>new Promise(resolve=>requests.push({args,resolve})) };
  const store=load(`src/store/modules/${kind}-storage.store.js`,kind==="experiment"?"ExperimentService":"ViewService",service);
  const context={commit:(name,value)=>store.mutations[name](store.state,value)};
  const first=store.actions["fetch"+plural](context,{page:1,pageSize:10,search:"alpha"});
  const second=store.actions["fetch"+plural](context,{page:1,pageSize:10,search:"beta"});
  assert.equal(requests.length,2);
  assert.equal(requests[0].args.search,"alpha");
  requests[1].resolve({count:1,results:[{[id]:2,name:"beta"}]}); await second;
  requests[0].resolve({count:20,results:[{[id]:1,name:"alpha"}]}); await first;
  const getters={}; for(const [key,value] of Object.entries(store.getters)) getters[key]=value(store.state,getters);
  const getter=kind==="view"?"getViewsPage":"getExperiments";
  const beta={page:1,pageSize:10,search:"beta"},alpha={...beta,search:"alpha"};
  assert.deepEqual(getters[getter](beta).map(row=>row[id]),[2]);
  assert.equal(getters["get"+plural+"Pagination"](beta).total,1);
  assert.equal(getters["get"+plural+"Pagination"](alpha).total,20);
  assert.equal(getters[getter]({...beta,search:""}),null);
 }
})().catch(error=>{console.error(error);process.exitCode=1;});
'''
    result = subprocess.run(["node", "-e", script], cwd=Path(__file__).resolve().parents[1], text=True, capture_output=True)
    assert result.returncode == 0, result.stderr


def test_changing_search_resets_pagination_and_hidden_selections():
    script = r'''
const fs=require("fs"),assert=require("assert"),Vue=require("vue");
const source=fs.readFileSync("src/mixins/collection-search.js","utf8").replace("export default","return");
const mixin=new Function(source)(),requests=[];
(async()=>{
 const vm=new Vue({mixins:[mixin],collectionSelectionKey:"selectedIds",collectionFetchAction:"fetch",
  data:()=>({selectedIds:{}}),beforeCreate(){this.$store={dispatch:async(action,query)=>{requests.push(query);}};}});
 vm.page=3; await Vue.nextTick();
 vm.selectedIds={hidden:true}; vm.search=" Helium "; await Vue.nextTick();
 assert.equal(vm.page,1); assert.deepEqual(vm.selectedIds,{});
 assert.deepEqual(requests[requests.length-1],{page:1,pageSize:10,search:"Helium"});
 vm.search=""; await Vue.nextTick();
 assert.equal(requests[requests.length-1].search,"");
 vm.$destroy();
})().catch(error=>{console.error(error);process.exitCode=1;});
'''
    result = subprocess.run(["node", "-e", script], cwd=Path(__file__).resolve().parents[1], text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
