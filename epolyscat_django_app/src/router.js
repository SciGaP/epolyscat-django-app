import Vue from 'vue';
import Router from 'vue-router';
import CreateExperiment from "@/components/Pages/CreateExperiment";
import Experiments from "@/components/Pages/Experiments";
import Runs from "@/components/Pages/Runs";
import ResubmitRun from "@/components/Pages/ResubmitRun";
import Views from "@/components/Pages/Views";
import Run from "@/components/Pages/Run";
import Viewable from "@/components/Pages/Viewable";

// Containers
//const Home = () => import('@/components/home');

const Login = () => import('@/components/Pages/Login');
const SignUp = () => import('@/components/Pages/SignUp');
const Home = () => import('@/components/Pages/Home');
const Documentation = () => import('@/components/Pages/Documentation');
const CreateRun = () => import('@/components/Pages/CreateRun')
const WorkflowRun = () => import('@/components/Pages/WorkflowRun')
const ViewRun = () => import('@/components/Pages/ViewRun')


Vue.use(Router);

export default new Router({
    mode: 'hash', // https://router.vuejs.org/api/#mode
    linkActiveClass: 'active',
    scrollBehavior: () => ({y: 0}),
    base: "/epolyscat_django_app/home/",
    routes: [
        {
            path: '/',
            name: 'home',
            component: Home,
        },
        {
            path: '/documentation',
            name: 'Documentation',
            component: Documentation,
        },
        {
            path: '/signup',
            name: 'SignUp',
            component: SignUp,
        },
        {
            path: '/login',
            name: 'Login',
            component: Login,
        },
        {
            path: '/experiments',
            name: 'Experiments',
            component: Experiments,
        },
        {
            path: '/views',
            name: 'Views',
            component: Views,
        },
        {
            path: '/runs',
            name: 'Runs',
            component: Runs,
        },
        {
            path: '/runs/new',
            name: 'CreateRun',
            component: CreateRun,
        },
        {
            path: '/runs/new/workflow',
            name: 'WorkflowRun',
            component: WorkflowRun,
        },
        {
            path: '/runs/new/legacy',
            name: 'LegacyCreateRun',
            component: Run,
        },

        {
            path: '/runs/:runId',
            name: 'Run',
            component: ViewRun,
        },
        {
            path: '/runs/:runId/legacy',
            name: 'LegacyRun',
            component: Run,
        },
        {
            path: '/runs/:runId/viewables/:filename',
            name: 'Viewable',
            component: Viewable,
        },
        {
            path: '/create-experiment',
            name: 'Create Experiment',
            component: CreateExperiment,
        },
        {
            path: '/create-run',
            name: 'CreateRunAlias',
            component: CreateRun
        },
        {
            path: '/views/:viewId',
            name: 'ViewRuns',
            component: Runs,
        },
        {
            path: '/tutorials',
            name: 'TutorialRuns',
            component: Runs
        },
        {
            path: '/resubmit-run',
            name: 'ResubmitRun',
            component: ResubmitRun
        }
    ],
});
