1. [Sign up](https://lamin.ai/signup) for a free account (see more [info](https://lamin.ai/docs/setup)) and copy the API key.
2. Log in on the command line:
   ```shell
   lamin login <email> --key <API-key>
   ```

```{note}

An account is free & [signing up](https://lamin.ai/signup) takes 1 min.

Through signing up, Lamin does _not_ store or see any of your data, but only _basic_ metadata about you (email address, etc.).

If you register a LaminDB instance to LaminHub, we store metadata about the involved infrastructure (database server, storage locations, etc.). But we don't store any secrets when you call `lamin init` and can't see your data or metadata. If you want to connect your instance to LaminHub for exploration with read-write access: please [reach out](https://lamin.ai/contact).

For more, see [doc](inv:docs#access), [the source code](https://github.com/laminlabs/lamindb-setup), or the [privacy policy](https://lamin.ai/legal/privacy-policy).

```

On the command line, you can log in with either email or handle:

```
lamin login testuser1@lamin.ai
lamin login testuser1
```

If you don't have a cached API-key in your environment, you need to copy it from your lamin.ai account and pass it:

```
lamin login <email> --key <API-key>
```

Log out:

```
lamin login --logout
```
