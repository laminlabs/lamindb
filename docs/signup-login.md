Why do I have to sign up? → Data flow needs a user identity to answer questions like: Who modified which data when? Who shares this with me?

An account is free & [signing up](https://lamin.ai/signup) takes 1 min.

```{note}

Lamin does _not_ store or see any of your data, but only _basic_ metadata about you (email address, etc.) & your LaminDB instances (S3 bucket names, etc.): see the [open-source client](https://github.com/laminlabs/lamindb-setup) & the [privacy policy](https://lamin.ai/legal/privacy-policy).

```

On the command line, you can log in with either email or handle:

```
lamin login testuser1@lamin.ai
lamin login testuser1
```

If you don't have a cached key in your environment, you need to pass it:

```
lamin login <email> --key <key>
```
