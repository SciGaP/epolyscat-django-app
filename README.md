

# ePolyscat Django App

## Getting Started

Follow the instructions for installing the [Airavata Django Portal](https://github.com/apache/airavata-django-portal).

With the Airavata Django Portal virtual environment activated, clone this repo and install it into the portal's virtual environment:

```bash
git clone --recursive https://github.com/SciGaP/epolyscat-django-app.git
cd epolyscat-django-app
pip install -e .
````

Start (or restart) the Django Portal server.

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