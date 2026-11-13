

# ePolyscat Django App

## Getting Started

Follow the instructions for installing the [Airavata Django Portal](https://github.com/apache/airavata-django-portal).

With the Airavata Django Portal virtual environment activated, clone this repo and install it into the portal's virtual environment:

```bash
git clone --recursive https://github.com/SciGaP/epolyscat-django-app.git
cd epolyscat-django-app
export NODE_OPTIONS=--openssl-legacy-provider
pip install -e .
````

Start (or restart) the Django Portal server.

### Application discovery

The app fetches the current user's accessible Airavata application modules, interfaces, and deployments before showing supported run choices. Gaussian16, OpenMolcas, and ePolyScat are matched by normalized module name; workflow stages and input mappings remain local application contracts. Utilities are read from the shared ePolyScat interface's `Application_Utility` editor options and checked against those contracts.

Optional `EPOLYSCAT_APPLICATION_ID`, `GAUSSIAN16_APPLICATION_ID`, and `OPENMOLCAS_APPLICATION_ID` entries in Django's `EPOLYSCAT` setting can select a particular registration. Configured IDs must still exist and be accessible. Missing or ambiguous registrations are unavailable, and discovery failures show a retry message rather than falling back to fixed IDs. Resource availability is checked separately for the selected allocation.

---

## Frontend Development

To build or develop the Vue.js frontend code, follow these steps:

1. Install **Node 16**
2. Install **Yarn**
3. In `./epolyscat_django_app/`, run:

```bash
cd ./epolyscat_django_app
yarn
yarn build
```

Or, to start the dev server for local development:

```bash
cd ./epolyscat_django_app
yarn
yarn serve
```

---

## Creating DB Migrations

```bash
django-admin makemigrations --pythonpath . --settings epolyscat_django_app.tests.settings epolyscat_django_app
```
