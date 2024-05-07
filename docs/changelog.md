# Changelog

:::{note}

If using LaminHub, please use the latest version of lamindb.

LaminDB implements "migration-based versioning". When upgrading your LaminDB installation to a new `minor` version in `major.minor.patch`, you also migrate your database by calling `lamin migrate deploy`.

:::

```{eval-rst}
.. role:: small
```

## 0.71

### 0.71.1 {small}`2024-05-07`

- 🚸  Upload source code of scripts upon `ln.finish()` and no longer upon `ln.track()` [PR](https://github.com/laminlabs/lamindb/pull/1624) [@falexwolf](https://github.com/falexwolf)
- 🎨 Make `features.add_feature_set` public [PR](https://github.com/laminlabs/lamindb/pull/1626) [@sunnyosun](https://github.com/sunnyosun)
- 🎨 Use the same `uid` for the same feature set in transfer [PR](https://github.com/laminlabs/lamindb/pull/1621) [@sunnyosun](https://github.com/sunnyosun)
- 🎨 Upon upload switch to virtual key [PR](https://github.com/laminlabs/lamindb/pull/1622) [@falexwolf](https://github.com/falexwolf)
- ⚡️ Zarr and cache improvements [PR](https://github.com/laminlabs/lamindb/pull/1620) [@Koncopd](https://github.com/Koncopd)
- ♻️ Extend valid suffixes to composite suffixes [PR](https://github.com/laminlabs/lamindb/pull/1619) [@falexwolf](https://github.com/falexwolf)
- 🔥 Remove little-used `artifact.view_tree()` [PR](https://github.com/laminlabs/lamindb/pull/1627) [@falexwolf](https://github.com/falexwolf)

### 0.71.0 {small}`2024-05-01`

- ✨ Manage multiple storage locations with integrity [PR](https://github.com/laminlabs/lamindb/pull/1611) [@falexwolf](https://github.com/falexwolf)
- 🚚 Add an `instance_uid` field to `Storage` | [374](https://github.com/laminlabs/lnschema-core/pull/374) [falexwolf](https://github.com/falexwolf)
- 🚸 Proper progress bars for upload and download [PR](https://github.com/laminlabs/lamindb/pull/1610) [@Koncopd](https://github.com/Koncopd)
- 🚸 Make save return self [PR](https://github.com/laminlabs/lamindb/pull/1606) [@falexwolf](https://github.com/falexwolf)

## 0.70

### 0.70.4 {small}`2024-04-24`

- ✨ Allow passing path to `.from_anndata` [PR](https://github.com/laminlabs/lamindb/pull/1600) [@sunnyosun](https://github.com/sunnyosun)
- 🚸 In `.setup.delete()`, check for data deletion & delete from hub [PR](https://github.com/laminlabs/lamindb/pull/1595) [@falexwolf](https://github.com/falexwolf)
- ⚡️ Speed up `latest_version` [PR](https://github.com/laminlabs/lamindb/pull/1594) [@falexwolf](https://github.com/falexwolf)
- 🚸 Better user feedback on folder-like artifacts [PR](https://github.com/laminlabs/lamindb/pull/1589) [@falexwolf](https://github.com/falexwolf)

### 0.70.3 {small}`2024-04-22`

- 🚸 Update metadata like description upon re-running [PR](https://github.com/laminlabs/lamindb/pull/1588) [@falexwolf](https://github.com/falexwolf)
- 🐛 Fix detection of AnnData in zarr and h5ad, refactor directory upload [PR](https://github.com/laminlabs/lamindb/pull/1587) [@Koncopd](https://github.com/Koncopd)
- 🚸 Raise error if transforms of type notebook or script are passed manually [PR](https://github.com/laminlabs/lamindb/pull/1584) [@falexwolf](https://github.com/falexwolf)

### 0.70.2 {small}`2024-04-19`

- ♻️ In Vitessce integration, separate `VitessceConfig` from its referenced artifacts [PR](https://github.com/laminlabs/lamindb/pull/1582) [@falexwolf](https://github.com/falexwolf)
- 🚸 In `ln.finish()`, remove flag `i_saved_the_notebook` [PR](https://github.com/laminlabs/lamindb/pull/1581) [@falexwolf](https://github.com/falexwolf)

### 0.70.1 {small}`2024-04-18`

- 🐛 Fix `public_source` in inspect [PR](https://github.com/laminlabs/lamindb/pull/1578) [@sunnyosun](https://github.com/sunnyosun)

### 0.70.0 {small}`2024-04-17`

- 🚸 Update data source in case transform is re-run [PR](https://github.com/laminlabs/lamindb/pull/1571) [@falexwolf](https://github.com/falexwolf)
- 🚸 Enable to label transforms via `transform.ulabels` [PR](https://github.com/laminlabs/lnschema-core/pull/370) [@falexwolf](https://github.com/falexwolf)
- 🚚 Deprecate `stage()` in favor of `cache()` [PR](https://github.com/laminlabs/lamindb/pull/1572) [@falexwolf](https://github.com/falexwolf)

## 0.69

### 0.69.10 {small}`2024-04-12`

- ✨ Add `.obsm` and  `.layers` to `MappedCollection` and rename `label_keys` to `obs_keys` [PR](https://github.com/laminlabs/lamindb/pull/1562) [@Koncopd](https://github.com/Koncopd)
- 🚸 Eliminate kwargs [PR](https://github.com/laminlabs/lamindb/pull/1561) [@sunnyosun](https://github.com/sunnyosun)
- ✨ Introduce `Annotate.from_mudata` [PR](https://github.com/laminlabs/lamindb/pull/1554) [@sunnyosun](https://github.com/sunnyosun)

### 0.69.9 {small}`2024-04-08`

- 🐛 Fix clashes for multiple processes [PR](https://github.com/laminlabs/lamindb/pull/1553) [@falexwolf](https://github.com/falexwolf)

### 0.69.8 {small}`2024-04-04`

- ♻️ Use future annotations [PR](https://github.com/laminlabs/lamindb/pull/1549) [@Zethson](https://github.com/Zethson)

### 0.69.7 {small}`2024-04-03`

- ✨ Add ability to upload arbitrary files or folders from CLI [PR](https://github.com/laminlabs/lamindb/pull/1545) [@falexwolf](https://github.com/falexwolf)
- 🐛 Fix anndata backed mode incompatibility with scipy 1.13.0 f

### 0.69.6 {small}`2024-04-02`

- 🚑️ Temp fix region for non-hosted buckets [PR](https://github.com/laminlabs/lamindb/pull/1543) [@sunnyosun](https://github.com/sunnyosun)

### 0.69.5 {small}`2024-03-30`

- ♻️ Improve Annotate API [PR](https://github.com/laminlabs/lamindb/pull/1542) [PR](https://github.com/laminlabs/lamindb/pull/1539)  [@sunnyosun](https://github.com/sunnyosun) [@falexwolf](https://github.com/falexwolf)
- ✨ Introduce `Registry.get()` and `lamin get` (replaces `lamin stage`) [PR](https://github.com/laminlabs/lamindb/pull/1538) [@falexwolf](https://github.com/falexwolf)

### 0.69.4 {small}`2024-03-30`

- ♻️ Add Vitessce integration [PR](https://github.com/laminlabs/lamindb/pull/1532) [@falexwolf](https://github.com/falexwolf)
- ♻️ Refactor collections [PR](https://github.com/laminlabs/lamindb/pull/1531) [@falexwolf](https://github.com/falexwolf)

### 0.69.3 {small}`2024-03-28`

- ✨ Introduce annotation flow via `Annotate.from_df` and `Annotate.from_anndata` [PR 1](https://github.com/laminlabs/lamindb/pull/1524) [2](https://github.com/laminlabs/lamindb/pull/1526) [3](https://github.com/laminlabs/lamindb/pull/1528) [@sunnyosun](https://github.com/sunnyosun)

### 0.69.2 {small}`2024-03-26`

- ✨ Stage collections [PR](https://github.com/laminlabs/lamindb/pull/1521) [@Koncopd](https://github.com/Koncopd)
- ✨ Improve functionality for folder-like artifacts [PR](https://github.com/laminlabs/lamindb/pull/1517) [@Koncopd](https://github.com/Koncopd)
- 📝 Improve the introduction page [PR](https://github.com/laminlabs/lamindb/pull/1510) [PR](https://github.com/laminlabs/lamindb/pull/1514) [@sunnyosun](https://github.com/sunnyosun) [@Zethson](https://github.com/Zethson)

### 0.69.1 {small}`2024-03-18`

✨ To try out, add `lamindb.validation` with the `Validator` class [PR](https://github.com/laminlabs/lamindb/pull/1508) [@sunnyosun](https://github.com/sunnyosun)

### 0.69.0 {small}`2024-03-17`

Main new features:

- ✨ Integrate lamindb with git [PR](https://github.com/laminlabs/lamindb/pull/1493) [PR](https://github.com/laminlabs/lamindb/pull/1497) [@falexwolf](https://github.com/falexwolf)
- ✨ Introduce `ln.finish()`, track run finish times as `run.finished_at`, rename `run.run_at` to `run.started_at`, upload notebooks during `ln.finish()` [PR](https://github.com/laminlabs/lamindb/pull/1501) [@falexwolf](https://github.com/falexwolf)
- 🚸 Upload script source code and environment during `ln.track()` [PR](https://github.com/laminlabs/lamindb/pull/1499) [@falexwolf](https://github.com/falexwolf)

Other changes:

- ✨ Allow including simple related fields in `.df()` [PR](https://github.com/laminlabs/lamindb/pull/1495) [@falexwolf](https://github.com/falexwolf)
- 🚚 Move transform settings into settings [PR](https://github.com/laminlabs/lamindb/pull/1498) [@falexwolf](https://github.com/falexwolf)
- ✨ Add `latest_version` filter for `QuerySet` [PR](https://github.com/laminlabs/lamindb/pull/1489) [@falexwolf](https://github.com/falexwolf)
- 🚚 Rename `transform.short_name` to `transform.key` [PR](https://github.com/laminlabs/lamindb/pull/1500) [@falexwolf](https://github.com/falexwolf)
- 🚸 Return `storage_idx` in `MappedCollection` [PR](https://github.com/laminlabs/lamindb/pull/1504) [@Koncopd](https://github.com/Koncopd)
- ♻️ Add a JSON field to `Run` [PR](https://github.com/laminlabs/lamindb/pull/1505) [@falexwolf](https://github.com/falexwolf)

## 0.68

### 0.68.2 {small}`2024-03-11`

- 🚸 Move transform & run artifacts into cache before uploading [PR](https://github.com/laminlabs/lamindb/pull/1488) [@falexwolf](https://github.com/falexwolf)
- 🚸 More sensible transform types [PR](https://github.com/laminlabs/lamindb/pull/1486) [@falexwolf](https://github.com/falexwolf)
- 🚚 Rename `lnschema_lamin1` to `wetlab` [PR](https://github.com/laminlabs/lamindb/pull/1487) [@falexwolf](https://github.com/falexwolf)

### 0.68.1 {small}`2024-03-08`

- 🚸 You can now use `ln.connect()` to connect to a LaminDB instance [PR](https://github.com/laminlabs/lamindb/pull/1480) [@falexwolf](https://github.com/falexwolf)
- 🚸 You can no longer delete data from non-default storage locations, as these might be tracked in other instances [PR](https://github.com/laminlabs/lamindb/pull/1484) [@sunnyosun](https://github.com/sunnyosun)
- 🚸 Enable transferring data from local instances to remote instances [PR](https://github.com/laminlabs/lamindb/pull/1479) [@sunnyosun](https://github.com/sunnyosun)

### 0.68.0 {small}`2024-03-01`

🚸 Decouple features linking from Artifact construction [PR 1](https://github.com/laminlabs/lamindb/pull/1434) [2](https://github.com/laminlabs/lamindb/pull/1455) [3](https://github.com/laminlabs/lamindb/pull/1458) [@sunnyosun](https://github.com/sunnyosun).

```python
# default constructor for PathLike
artifact = ln.Artifact("mysc.h5ad", description="raw data")
# from_ constructors for other types
artifact = ln.Artifact.from_anndata(mysc_adata, description="raw data")  # no longer links features
artifact = artifact.save()

# high-level feature linking
artifact.features.add_from_anndata(var_field=bt.Gene.ensembl_gene_id)
artifact.features.add_from_df()

# low-level feature linking
meta = ln.Feature.from_values(mysc_adata.obs.columns, field="name")
genes = bt.Gene.from_values(mysc_adata.var.ensembl_gene_id, field="ensembl_gene_id")
artifact.features.add(genes, slot="obs")
artifact.features.add(genes, slot="var")

# labels linking (no change)
labels = ln.ULabel.from_values(adata.obs.donor, field=...)
ln.save(labels)
artifact.labels.add(labels)
```

<br>

- 🚸 Can now use `ln.track()` without `lamin track` [PR](https://github.com/laminlabs/lamindb/pull/1462) [@falexwolf](https://github.com/falexwolf)
- 🐛 `lamin stage` respects new URL design [PR](https://github.com/laminlabs/lamindb/pull/1467) [@falexwolf](https://github.com/falexwolf)
- 🚚 Rename `.dev` to `.core` [PR](https://github.com/laminlabs/lamindb/pull/1464) [@falexwolf](https://github.com/falexwolf)
- ♻️ Improved `MappedCollection` [PR](https://github.com/laminlabs/lamindb/pull/1460) [PR](https://github.com/laminlabs/lamindb/pull/1448) [@Koncopd](https://github.com/Koncopd)


## 0.67

### 0.67.3 {small}`2024-02-02`

- 🚸 Can now import `bionty` instead of `lnschema-bionty` [PR](https://github.com/laminlabs/lamindb/pull/1415) [@sunnyosun](https://github.com/sunnyosun)
- ♻️ Use Click for CLI [PR](https://github.com/laminlabs/lamindb/pull/1420) [@chaichontat](https://github.com/chaichontat)
- 🎨 Make `Collection.save()` ACID [PR](https://github.com/laminlabs/lamindb/pull/1410) [@falexwolf](https://github.com/falexwolf)
- ✨ Add `add_to_version_family` [PR](https://github.com/laminlabs/lamindb/pull/1408) [@sunnyosun](https://github.com/sunnyosun)
- 🐛 Transfer `collection.artifacts` [PR](https://github.com/laminlabs/lamindb/pull/1405) [@sunnyosun](https://github.com/sunnyosun)

### 0.67.2 {small}`2024-01-14`

- ✨ Enable staging notebooks & code using the CLI [PR](https://github.com/laminlabs/lamindb/pull/1403) [@falexwolf](https://github.com/falexwolf)

### 0.67.1 {small}`2024-01-12`

- 🐛 Fix idempotency of `collection.save()` [PR](https://github.com/laminlabs/lamindb/pull/1401) [@falexwolf](https://github.com/falexwolf)
- 🚸 Disallow bulk-delete for Artifact, Transform & Collection [PR](https://github.com/laminlabs/lamindb/pull/1398) [@falexwolf](https://github.com/falexwolf)
- 🚸 Init transform versions at 1 [PR](https://github.com/laminlabs/lamindb/pull/1397) [@falexwolf](https://github.com/falexwolf)
- ✨ Load json and html files [PR](https://github.com/laminlabs/lamindb/pull/1396) [@falexwolf](https://github.com/falexwolf)

### 0.67.0 {small}`2024-01-11`

- 🚚 Rename `.bionty` to `.public`, `.from_bionty` to `.from_public` [PR](https://github.com/laminlabs/lamindb/pull/1394) [@sunnyosun](https://github.com/sunnyosun)

## 0.66

### 0.66.1 {small}`2024-01-09`

- 🐛 Fix id matching in view_lineage [PR](https://github.com/laminlabs/lamindb/pull/1395) [@sunnyosun](https://github.com/sunnyosun)
- ♻️ Fix connection time outs [PR](https://github.com/laminlabs/lamindb-setup/pull/611) [@Koncopd](https://github.com/Koncopd)
- ♻️ Incorporate edge cases in `inner` and `outer` join in `Collection.mapped` [PR](https://github.com/laminlabs/lamindb/pull/1392) [@Koncopd](https://github.com/Koncopd)
- 🎨 Not create organism records when calling `.bionty()` [PR](https://github.com/laminlabs/lamindb/pull/1391) [@sunnyosun](https://github.com/sunnyosun)

### 0.66.0 {small}`2024-01-07`

- 🚸 Add anonymous access (now works without login) [PR1386](https://github.com/laminlabs/lamindb/pull/1386) [@falexwolf](https://github.com/falexwolf)
- 🎨 Introduce ordered collections and simplify `.mapped()` [PR1390](https://github.com/laminlabs/lamindb/pull/1390) [@falexwolf](https://github.com/falexwolf)
- 📝 Re-write quickstart [PR1387](https://github.com/laminlabs/lamindb/pull/1387) [@falexwolf](https://github.com/falexwolf)

## 0.65

### 0.65.1 {small}`2024-01-05`

- 🩹 Prepare a potential migration of the hub to Django [PR1385](https://github.com/laminlabs/lamindb/pull/1385) [@falexwolf](https://github.com/falexwolf)
- 🚸 Various improvements [PR1384](https://github.com/laminlabs/lamindb/pull/1384) [@falexwolf](https://github.com/falexwolf)
- 🩹 Track suffix of requirements.txt [PR1383](https://github.com/laminlabs/lamindb/pull/1383) [@falexwolf](https://github.com/falexwolf)
- ✨ Add outer join and categories caching to Collection.mapped [PR1380](https://github.com/laminlabs/lamindb/pull/1380) [@Koncopd](https://github.com/Koncopd)
- ♻️ Except memory error [PR1382](https://github.com/laminlabs/lamindb/pull/1382) [@falexwolf](https://github.com/falexwolf)

### 0.65.0 {small}`2024-01-02`

- 🚚 Rename `Dataset` to `Collection` [PR1377](https://github.com/laminlabs/lamindb/pull/1377) [@falexwolf](https://github.com/falexwolf)
- ✨ Track run environment [PR1368](https://github.com/laminlabs/lamindb/pull/1368) [@falexwolf](https://github.com/falexwolf)
- ✨ Allow transfer from private instances [PR1370](https://github.com/laminlabs/lamindb/pull/1370) [@falexwolf](https://github.com/falexwolf)
- 🚸 Speed up transfer and enable transfer parents [PR1371](https://github.com/laminlabs/lamindb/pull/1371) [@sunnyosun](https://github.com/sunnyosun)
- 🎨 Version based on `stem_uid` instead of `initial_version_id` and replace `__lamindb_uid_prefix__` with `__transform_stem_uid__` [PR1369](https://github.com/laminlabs/lamindb/pull/1369) [PR1375](https://github.com/laminlabs/lamindb/pull/1375) [PR1373](https://github.com/laminlabs/lamindb/pull/1373) [@bpenteado](https://github.com/bpenteado) [@falexwolf](https://github.com/falexwolf)
- 🎨 Name `.lndb` files by instance id [PR1372](https://github.com/laminlabs/lamindb/pull/1372) [@falexwolf](https://github.com/falexwolf)

## 0.64

### 0.64.2 {small}`2023-12-18`

- 🚚 Deprecate `view_flow` in favor of `view_lineage` [PR1367](https://github.com/laminlabs/lamindb/pull/1367) [@falexwolf](https://github.com/falexwolf)
- 🐛 Account for length zero label lists [PR1362](https://github.com/laminlabs/lamindb/pull/1362) [@falexwolf](https://github.com/falexwolf)
- 🐛 Fix mypy [PR1359](https://github.com/laminlabs/lamindb/pull/1359) [@Koncopd](https://github.com/Koncopd)
- 👷 Use ruff and clean up pre-commit [PR1358](https://github.com/laminlabs/lamindb/pull/1358) [@Zethson](https://github.com/Zethson)

### 0.64.1 {small}`2023-12-13`

- 🚸 Silence the httpx logger [PR1356](https://github.com/laminlabs/lamindb/pull/1356) [@falexwolf](https://github.com/falexwolf)
- 🐛 Fix hash inference for cloud non-s3 objects and use stored gcs credentials [PR1355](https://github.com/laminlabs/lamindb/pull/1355) [@Koncopd](https://github.com/Koncopd)

### 0.64.0 {small}`2023-12-13`

🏗️ Replace `File` with `Artifact` to model directories and directory-like array stores [PR1345](https://github.com/laminlabs/lamindb/pull/1345) [@falexwolf](https://github.com/falexwolf).

For artifacts that are simple files, rename `File` to `Artifact` along with the word `file` to `artifact` in registry fields. Otherwise, your code will behave the same.

## 0.63

### 0.63.5 {small}`2023-12-11`

- 🐛 Fix `lamin save` for scripts [PR1346](https://github.com/laminlabs/lamindb/pull/1346) [@sunnyosun](https://github.com/sunnyosun)

### 0.63.4 {small}`2023-12-08`

- 🚸 Parse variables upon track, do not access storage during load [PR1341](https://github.com/laminlabs/lamindb/pull/1341) [@falexwolf](https://github.com/falexwolf)
- 🚸 Raise error when tracking script fails [PR1340](https://github.com/laminlabs/lamindb/pull/1340) [@sunnyosun](https://github.com/sunnyosun)
- 🚸 Literal checks for `type` when creating a Feature [PR1339](https://github.com/laminlabs/lamindb/pull/1339) [@sunnyosun](https://github.com/sunnyosun)
- 🐛 Fix capitalization of feature names in describe [PR1337](https://github.com/laminlabs/lamindb/pull/1337) [@sunnyosun](https://github.com/sunnyosun)

### 0.63.3 {small}`2023-12-05`

- ✨ Do virtual inner join of variables in Dataset.mapped [PR1334](https://github.com/laminlabs/lamindb/pull/1334) [@Koncopd](https://github.com/Koncopd)
- 📝 Add a track notebooks how-to [PR1333](https://github.com/laminlabs/lamindb/pull/1333) [@sunnyosun](https://github.com/sunnyosun)
- 🚚 Move high-level bionty docs here [PR1332](https://github.com/laminlabs/lamindb/pull/1332) [@falexwolf](https://github.com/falexwolf)

### 0.63.2 {small}`2023-12-02`

- ✨ Add `return_field` to `standardize()` [PR1331](https://github.com/laminlabs/lamindb/pull/1331) [@sunnyosun](https://github.com/sunnyosun)
- ♻️ Extend `MappedDataset.get_label_weights` [PR1330](https://github.com/laminlabs/lamindb/pull/1330) [@Koncopd](https://github.com/Koncopd)

### 0.63.1 {small}`2023-11-28`

- 🚸 Enable schema packages that do not start with `lnschema_` [592](https://github.com/laminlabs/lamindb-setup/pull/592) [falexwolf](https://github.com/falexwolf)
- 🚸 Do not rely on naming the primary key id [PR1328](https://github.com/laminlabs/lamindb/pull/1328) [@falexwolf](https://github.com/falexwolf)
- ✨ Add context manager to `File.backed` and `Dataset.mapped` [PR1326](https://github.com/laminlabs/lamindb/pull/1326) [@Koncopd](https://github.com/Koncopd)

### 0.63.0 {small}`2023-11-24`

- ♻️ Centralize progressbar and convert to percentage [PR1322](https://github.com/laminlabs/lamindb/pull/1322) [@falexwolf](https://github.com/falexwolf)
- ✨ Added view_tree for keys [PR1321](https://github.com/laminlabs/lamindb/pull/1321) [@sunnyosun](https://github.com/sunnyosun)
- 🚚 Recode visibility to default 1, hidden 0, trash -1 [PR1319](https://github.com/laminlabs/lamindb/pull/1319) [@sunnyosun](https://github.com/sunnyosun)
- ✨ Track Python scripts [PR1318](https://github.com/laminlabs/lamindb/pull/1318) [@falexwolf](https://github.com/falexwolf)
- 📝 Modify Quickstart to use File instead of Dataset [PR1317](https://github.com/laminlabs/lamindb/pull/1317) [@sunnyosun](https://github.com/sunnyosun)
- 📝 Add an FAQ on storage keys & locations [PR1253](https://github.com/laminlabs/lamindb/pull/1253) [@Zethson](https://github.com/Zethson)

## 0.62

### 0.62.0 {small}`2023-11-20`

- ✨ Create `dataset.mapped()` for data loaders over `AnnData` [PR1307](https://github.com/laminlabs/lamindb/pull/1307) [@Koncopd](https://github.com/Koncopd)
- 🚸 Add `.df()` to `Registry` [PR1313](https://github.com/laminlabs/lamindb/pull/1313) [@falexwolf](https://github.com/falexwolf)
- 🐛 Fix passing organism to featureset.from_values [PR1312](https://github.com/laminlabs/lamindb/pull/1312) [@sunnyosun](https://github.com/sunnyosun)
- 🐛 Fix idempotency of transfer files [PR1311](https://github.com/laminlabs/lamindb/pull/1311) [@falexwolf](https://github.com/falexwolf)
- 🍱 Updated bionty sources [PR1305](https://github.com/laminlabs/lamindb/pull/1305) [@sunnyosun](https://github.com/sunnyosun)
- 📝 Add an annotation guide [PR1302](https://github.com/laminlabs/lamindb/pull/1302) [@sunnyosun](https://github.com/sunnyosun)

## 0.61

### 0.61.0 {small}`2023-11-15`

- 🚸 Prettify `view_tree`, check for ahead migrations [PR1299](https://github.com/laminlabs/lamindb/pull/1299) [@sunnyosun](https://github.com/sunnyosun)
- 🚸 Ensure saved notebooks aren't rerun with the same version [PR1298](https://github.com/laminlabs/lamindb/pull/1298) [@falexwolf](https://github.com/falexwolf)
- 🔥 Remove `Modality` registry [PR1297](https://github.com/laminlabs/lamindb/pull/1297) [@sunnyosun](https://github.com/sunnyosun)
- ⬆️ More flexible versions for anndata and s3fs [PR1295](https://github.com/laminlabs/lamindb/pull/1295) [@Koncopd](https://github.com/Koncopd)
- ♻️ Refactor `Dataset` construction from directory [PR1300](https://github.com/laminlabs/lamindb/pull/1300) [@falexwolf](https://github.com/falexwolf)

## 0.60

### 0.60.3 {small}`2023-11-13`

- 🚑 Fix visibility choices, updated intro [PR1294](https://github.com/laminlabs/lamindb/pull/1294) [@sunnyosun](https://github.com/sunnyosun)

### 0.60.1 {small}`2023-11-10`

- 📝 Clearer badge titles [PR1291](https://github.com/laminlabs/lamindb/pull/1291) [@falexwolf](https://github.com/falexwolf)
- ⏪ Fix bionty versionsdir [PR1290](https://github.com/laminlabs/lamindb/pull/1290) [@sunnyosun](https://github.com/sunnyosun)

### 0.60.0 {small}`2023-11-09`

- 🚸 Enable to pass metadata to dataset [PR1289](https://github.com/laminlabs/lamindb/pull/1289) [@falexwolf](https://github.com/falexwolf)
- 🚚 Move `view_tree` to `UPath` [PR1288](https://github.com/laminlabs/lamindb/pull/1288) [@falexwolf](https://github.com/falexwolf)
- ♻️ Move `UPath` to root-level API [PR1286](https://github.com/laminlabs/lamindb/pull/1286) [@falexwolf](https://github.com/falexwolf)
- 🚸 Enable to add labels without a feature [PR1285](https://github.com/laminlabs/lamindb/pull/1285) [@sunnyosun](https://github.com/sunnyosun)

## 0.59

### 0.59.3 {small}`2023-11-04`

- 💄 Show transform uid in view_lineage [PR1280](https://github.com/laminlabs/lamindb/pull/1280) [@sunnyosun](https://github.com/sunnyosun)
- 🚸 Improve search [PR1279](https://github.com/laminlabs/lamindb/pull/1279) [@sunnyosun](https://github.com/sunnyosun)

### 0.59.2 {small}`2023-11-04`

- 🐛 Backed raw subset fix [PR1278](https://github.com/laminlabs/lamindb/pull/1278) [@Koncopd](https://github.com/Koncopd)
- 🚸 Improve connection stability [PR1277](https://github.com/laminlabs/lamindb/pull/1277) [@falexwolf](https://github.com/falexwolf)
- 🚸 Improve search [PR1275](https://github.com/laminlabs/lamindb/pull/1275) [@sunnyosun](https://github.com/sunnyosun)

### 0.59.1 {small}`2023-11-03`

- ♻️ Move CLI from `lamindb-setup` here [PR1264](https://github.com/laminlabs/lamindb/pull/1264) [@falexwolf](https://github.com/falexwolf)
- 🐛 `file.replace()` works with virtual keys [PR1265](https://github.com/laminlabs/lamindb/pull/1265) [@Koncopd](https://github.com/Koncopd)

### 0.59.0 {small}`2023-10-27`

- 🚸 Introduce virtual keys [PR1262](https://github.com/laminlabs/lamindb/pull/1262) [@falexwolf](https://github.com/falexwolf)
- ♻️ Refactor visibility [PR1261](https://github.com/laminlabs/lamindb/pull/1261) [@falexwolf](https://github.com/falexwolf)
- ♻️ Turn `using` into a `Registry` method [PR1259](https://github.com/laminlabs/lamindb/pull/1259) [@falexwolf](https://github.com/falexwolf)

## 0.58

### 0.58.2 {small}`2023-10-26`

- 👷 Enable to read db from env variable [PR1258](https://github.com/laminlabs/lamindb/pull/1258) [@falexwolf](https://github.com/falexwolf)

### 0.58.1 {small}`2023-10-26`

- 🔒️ Re-enable read-only access for public instances [PR1257](https://github.com/laminlabs/lamindb/pull/1257) [@falexwolf](https://github.com/falexwolf)
- ✨ Introduce `visibility` to `File` and `Dataset` [PR1254](https://github.com/laminlabs/lamindb/pull/1254) [@sunnyosun](https://github.com/sunnyosun)
- 🔒️ Enable to pass DB URLs during load [PR1256](https://github.com/laminlabs/lamindb/pull/1256) [@falexwolf](https://github.com/falexwolf)
- ✨ Enable passing filter expressions to lookup and search [PR1255](https://github.com/laminlabs/lamindb/pull/1255) [@sunnyosun](https://github.com/sunnyosun)
- ♻️ Use `file.uid` instead of `file.id` in `dataset.load()` [PR1250](https://github.com/laminlabs/lamindb/pull/1250) [@falexwolf](https://github.com/falexwolf)

### 0.58.0 {small}`2023-10-24`

- ⚡️ Speed up dataset creation in presence of millions of feature links [PR1248](https://github.com/laminlabs/lamindb/pull/1248) [@falexwolf](https://github.com/falexwolf)
- 📝 Add a faq for custom validator [PR1249](https://github.com/laminlabs/lamindb/pull/1249) [@sunnyosun](https://github.com/sunnyosun)
- ✨ Add `parents` to `Organism` [PR1245](https://github.com/laminlabs/lamindb/pull/1245) [@sunnyosun](https://github.com/sunnyosun)

## 0.57

### 0.57.2 {small}`2023-10-23`

- ➕ Allow more anndata versions [PR1244](https://github.com/laminlabs/lamindb/pull/1244) [@Koncopd](https://github.com/Koncopd)

### 0.57.1 {small}`2023-10-20`

- 🎨 Display uid instead of id in search result [PR1239](https://github.com/laminlabs/lamindb/pull/1239) [@sunnyosun](https://github.com/sunnyosun)

### 0.57.0 {small}`2023-10-19`

- 🚚 Rename `Species` to `Organism` [PR1237](https://github.com/laminlabs/lamindb/pull/1237) [@sunnyosun](https://github.com/sunnyosun)
- 🔥 Remove email field from `User` registry [PR1238](https://github.com/laminlabs/lamindb/pull/1238) [@falexwolf](https://github.com/falexwolf)

## 0.56

### 0.56.2 {small}`2023-10-19`

- 🚚 Migrate from password to apikey [PR1233](https://github.com/laminlabs/lamindb/pull/1233) [@falexwolf](https://github.com/falexwolf)

### 0.56.1 {small}`2023-10-18`

- 🐛 Fix transfer of features [PR1232](https://github.com/laminlabs/lamindb/pull/1232) [@sunnyosun](https://github.com/sunnyosun)

### 0.56.0 {small}`2023-10-17`

The previous universal primary keys are now stored us unique universal ids in `.uid` and `.id` stores integer primary keys.

This leads to significantly higher storage efficiency of indexes & constraints in link tables with >10 million rows.

A dialogue guides you through the migration process.

- 🚚 Migrate to integer primary keys [PR1229](https://github.com/laminlabs/lamindb/pull/1229) [@falexwolf](https://github.com/falexwolf)
- ♻️ More flexible cache management [PR1230](https://github.com/laminlabs/lamindb/pull/1230) [@Koncopd](https://github.com/Koncopd)

## 0.55

### 0.55.2 {small}`2023-10-10`

- ♻️ Keep raw in pbmc68k_reduced and human_immune_cells [PR1228](https://github.com/laminlabs/lamindb/pull/1228) [@Koncopd](https://github.com/Koncopd)
- ⬆️ Make backed compatible with anndata==0.10.1 [PR1227](https://github.com/laminlabs/lamindb/pull/1227) [@Koncopd](https://github.com/Koncopd)
- ✨ Also transfer features [PR1226](https://github.com/laminlabs/lamindb/pull/1226) [@falexwolf](https://github.com/falexwolf)

### 0.55.1 {small}`2023-10-4`

- ✨ Dataset from paths aka storage [PR1222](https://github.com/laminlabs/lamindb/pull/1222) [@falexwolf](https://github.com/falexwolf)
- 🚸 Do not search for alternative names if id is passed [PR1225](https://github.com/laminlabs/lamindb/pull/1225) [@falexwolf](https://github.com/falexwolf)
- 📝 Added lamindata badge [PR1223](https://github.com/laminlabs/lamindb/pull/1223) [@sunnyosun](https://github.com/sunnyosun)
- ♻️ Refactor data flow tracking for Dataset [PR1224](https://github.com/laminlabs/lamindb/pull/1224) [@falexwolf](https://github.com/falexwolf)

### 0.55.0 {small}`2023-10-03`

- ✨ Enable to transfer data between instances [PR1216](https://github.com/laminlabs/lamindb/pull/1216) [@falexwolf](https://github.com/falexwolf)
- ✨ Enable saving notebook files [PR1213](https://github.com/laminlabs/lamindb/pull/1213) [@falexwolf](https://github.com/falexwolf)
- 💄 Display id in view_parents [PR1215](https://github.com/laminlabs/lamindb/pull/1215) [@sunnyosun](https://github.com/sunnyosun)

## 0.54

### 0.54.4 {small}`2023-09-29`

- ✨ Add labels transfer through `labels.add_from` [PR1211](https://github.com/laminlabs/lamindb/pull/1211) [@falexwolf](https://github.com/falexwolf)
- ♻️ Perform outer join during load [PR1212](https://github.com/laminlabs/lamindb/pull/1212) [@falexwolf](https://github.com/falexwolf)
- 🎨 Show both file and dataset in view_lineage [PR1209](https://github.com/laminlabs/lamindb/pull/1209) [@sunnyosun](https://github.com/sunnyosun)

### 0.54.3 {small}`2023-09-27`

- 💄 Fix double arrows for file and dataset [PR1206](https://github.com/laminlabs/lamindb/pull/1206) [@sunnyosun](https://github.com/sunnyosun)
- 🎨 Add a members property to FeatureSet [PR1205](https://github.com/laminlabs/lamindb/pull/1205) [@falexwolf](https://github.com/falexwolf)

### 0.54.2 {small}`2023-09-26`

- 🎨 Better implementation of `dataset.view_lineage` [PR1194](https://github.com/laminlabs/lamindb/pull/1194) [@sunnyosun](https://github.com/sunnyosun)

### 0.54.1 {small}`2023-09-21`

- 🔇 Mute repetitive warning for saving parents [PR1190](https://github.com/laminlabs/lamindb/pull/1190) [@sunnyosun](https://github.com/sunnyosun)
- 🔇 Fix warning logging upon save [PR1185](https://github.com/laminlabs/lamindb/pull/1185) [@sunnyosun](https://github.com/sunnyosun)
- 🩹 Better suffix handling [PR1181](https://github.com/laminlabs/lamindb/pull/1181) [@falexwolf](https://github.com/falexwolf)

### 0.54.0 {small}`2023-09-17`

- ✏️ Fix views for interactive terminal [PR1178](https://github.com/laminlabs/lamindb/pull/1178) [@sunnyosun](https://github.com/sunnyosun)
- 🚸 Improve error message for conflicting storage keys [PR1174](https://github.com/laminlabs/lamindb/pull/1174) [@falexwolf](https://github.com/falexwolf)
- ♻️ Add string levels [PR1173](https://github.com/laminlabs/lamindb/pull/1173) [@falexwolf](https://github.com/falexwolf)
- ✨ Add `reference` and `reference_type` to `Dataset`, `ULabel` [PR1171](https://github.com/laminlabs/lamindb/pull/1171) [@sunnyosun](https://github.com/sunnyosun)
- 🔒 Add client fallbacks to lamindb-setup [PR1172](https://github.com/laminlabs/lamindb/pull/1172) [@falexwolf](https://github.com/falexwolf)
- 🐛 Fix `view_parents()` when there are no parents [PR1167](https://github.com/laminlabs/lamindb/pull/1167) [@sunnyosun](https://github.com/sunnyosun)
- 🔊 Turn on hint for tutorials [PR1164](https://github.com/laminlabs/lamindb/pull/1164) [@sunnyosun](https://github.com/sunnyosun)
- ⚡ Speed up instance loading time & refactor lamindb-setup [PR1163](https://github.com/laminlabs/lamindb/pull/1163) [@falexwolf](https://github.com/falexwolf)
- ⬆️ Update fsspec to 2023.9.0 [PR1162](https://github.com/laminlabs/lamindb/pull/1162) [@Koncopd](https://github.com/Koncopd)
- 🎨 More validation for `File` creation, improve zarr backend [PR1160](https://github.com/laminlabs/lamindb/pull/1160) [@Koncopd](https://github.com/Koncopd)
- 🔥 Remove `add_labels()` & `get_labels()` [PR1155](https://github.com/laminlabs/lamindb/pull/1155) [@falexwolf](https://github.com/falexwolf)
- 🎨 Increased search threshold to 90 [PR1156](https://github.com/laminlabs/lamindb/pull/1156) [@sunnyosun](https://github.com/sunnyosun)

## 0.53 {small}`2023-09-09`

- ♻️ Refactor hub client [PR1150](https://github.com/laminlabs/lamindb/pull/1150) [@falexwolf](https://github.com/falexwolf)
- ✨ Add `return_field` to `Lookup` [PR1153](https://github.com/laminlabs/lamindb/pull/1153) [@sunnyosun](https://github.com/sunnyosun)
- 🔇 Silence download bar for <5MB bionty files [PR1151](https://github.com/laminlabs/lamindb/pull/1151) [@sunnyosun](https://github.com/sunnyosun)
- 🔇 Set default verbosity to warning [PR1149](https://github.com/laminlabs/lamindb/pull/1149) [@sunnyosun](https://github.com/sunnyosun)
- 🚚 Add `Data.labels` accessor, rename `Label` to `ULabel` [PR1148](https://github.com/laminlabs/lamindb/pull/1148) [@sunnyosun](https://github.com/sunnyosun)
- ✨ Allow create records from non-default bionty sources [PR1147](https://github.com/laminlabs/lamindb/pull/1147) [@sunnyosun](https://github.com/sunnyosun)
- 🔥 Removed `get_by_feature` [PR1145](https://github.com/laminlabs/lamindb/pull/1145) [@sunnyosun](https://github.com/sunnyosun)
- 🚚 Replace `get_feature_set` with `.feature_sets.__getitem__` [PR1141](https://github.com/laminlabs/lamindb/pull/1141) [@sunnyosun](https://github.com/sunnyosun)

## 0.52

### 0.52.2 {small}`2023-09-05`

- 🎨 Add labels to describe [PR1139](https://github.com/laminlabs/lamindb/pull/1139) [@sunnyosun](https://github.com/sunnyosun)
- 🎨 Set default field for from_values [PR1137](https://github.com/laminlabs/lamindb/pull/1137) [@sunnyosun](https://github.com/sunnyosun)
- 📝 More comprehensive landing page [PR1136](https://github.com/laminlabs/lamindb/pull/1136) [@falexwolf](https://github.com/falexwolf)
- 🎨 Switch to WRatio and default to limit=20 for search [PR1135](https://github.com/laminlabs/lamindb/pull/1135) [@sunnyosun](https://github.com/sunnyosun)
- 🚚 Rename `File.tree` to `File.view_tree` [PR1133](https://github.com/laminlabs/lamindb/pull/1133) [@falexwolf](https://github.com/falexwolf)
- ✨ Enable dataset versioning [PR1132](https://github.com/laminlabs/lamindb/pull/1132) [@falexwolf](https://github.com/falexwolf)

### 0.52.1 {small}`2023-09-02`

- ♻️ Refactor features display [PR1129](https://github.com/laminlabs/lamindb/pull/1129) [@falexwolf](https://github.com/falexwolf)
- ♻️ Refactor feature set linking [PR1128](https://github.com/laminlabs/lamindb/pull/1128) [@falexwolf](https://github.com/falexwolf)
- ♻️ Add display of numerical features [PR1127](https://github.com/laminlabs/lamindb/pull/1127) [@falexwolf](https://github.com/falexwolf)
- 📝 Rewrite landing page [PR1126](https://github.com/laminlabs/lamindb/pull/1126) [@falexwolf](https://github.com/falexwolf)
- 🚸 Validate `Run` records [PR1123](https://github.com/laminlabs/lamindb/pull/1123) [@falexwolf](https://github.com/falexwolf)
- 🔊 Improved logging message for inspect [PR1124](https://github.com/laminlabs/lamindb/pull/1124) [@sunnyosun](https://github.com/sunnyosun)

### 0.52.0 {small}`2023-09-01`

- ♻️ Type features [PR1121](https://github.com/laminlabs/lamindb/pull/1121) [@falexwolf](https://github.com/falexwolf)
- 🚸 Enable to pass modality to `from_df` and `from_anndata` [PR1116](https://github.com/laminlabs/lamindb/pull/1116) [@falexwolf](https://github.com/falexwolf)
- ⬆️ Cache instance id in settings files [PR1120](https://github.com/laminlabs/lamindb/pull/1120) [@bpenteado](https://github.com/bpenteado)
- 🚑️ Fix map_synonyms bug [PR1119](https://github.com/laminlabs/lamindb/pull/1119) [@sunnyosun](https://github.com/sunnyosun)
- 🐛 Fix `target_label` in `view_parents` [PR1118](https://github.com/laminlabs/lamindb/pull/1118) [@sunnyosun](https://github.com/sunnyosun)
- ♻️ Simplify `.standardize()` implementation [PR1117](https://github.com/laminlabs/lamindb/pull/1117) [@sunnyosun](https://github.com/sunnyosun)
- 🐛 Fix `lookup()` for `QuerySet` [PR1115](https://github.com/laminlabs/lamindb/pull/1115) [@sunnyosun](https://github.com/sunnyosun)
- 🩹 Replace deprecated `is_categorical_dtype` [PR1114](https://github.com/laminlabs/lamindb/pull/1114) [@falexwolf](https://github.com/falexwolf)

## 0.51

### 0.51.3 {small}`2023-08-30`

- 🚸 Better treatment of `File.suffix` [PR1113](https://github.com/laminlabs/lamindb/pull/1113) [@falexwolf](https://github.com/falexwolf)

### 0.51.2 {small}`2023-08-28`

- ♻️ Enable typed feature sets for dataframes [PR1110](https://github.com/laminlabs/lamindb/pull/1110) [@falexwolf](https://github.com/falexwolf)
- ♻️ Refactor to use `UPath` everywhere [PR1102](https://github.com/laminlabs/lamindb/pull/1102) [@Koncopd](https://github.com/Koncopd)
- ⚡️ Delete cached files for local instances on File.save() [PR1109](https://github.com/laminlabs/lamindb/pull/1109) [@Koncopd](https://github.com/Koncopd)
- 🚚 Rename data lineage to data flow [PR1108](https://github.com/laminlabs/lamindb/pull/1108) [@falexwolf](https://github.com/falexwolf)

### 0.51.1 {small}`2023-08-28`

- 📝 Added `validate` guide [PR1104](https://github.com/laminlabs/lamindb/pull/1104) [@sunnyosun](https://github.com/sunnyosun)
- 🐛 Fix feature set updates & `describe()` if no features are present [PR1098](https://github.com/laminlabs/lamindb/pull/1098) [@falexwolf](https://github.com/falexwolf)
- 🎨 Make `species` not nullable in `Gene`, `Protein`, `CellType` [PR1095](https://github.com/laminlabs/lamindb/pull/1095) [@sunnyosun](https://github.com/sunnyosun)

### 0.51.0 {small}`2023-08-24`

- ♻️ Make `suffix` non-nullable and root unique [PR1090](https://github.com/laminlabs/lamindb/pull/1090) [@falexwolf](https://github.com/falexwolf)
- ♻️ Add `view_lineage` to `Dataset` [PR1089](https://github.com/laminlabs/lamindb/pull/1089) [@falexwolf](https://github.com/falexwolf)
- 🚸 Disable creating non-validated records with `from_values` [PR1086](https://github.com/laminlabs/lamindb/pull/1086) [@sunnyosun](https://github.com/sunnyosun)
- ♻️ Enable passing AnnData-like `File` object to `Dataset.from_anndata` [PR1085](https://github.com/laminlabs/lamindb/pull/1085) [@falexwolf](https://github.com/falexwolf)
- 🔊 More consistent logging of inspect [PR1081](https://github.com/laminlabs/lamindb/pull/1081) [@sunnyosun](https://github.com/sunnyosun)
- 🔥 Remove `inherit_relations` [PR1079](https://github.com/laminlabs/lamindb/pull/1079) [@sunnyosun](https://github.com/sunnyosun)
- ♻️ Use `CanValidate` & `HasParents` [PR1078](https://github.com/laminlabs/lamindb/pull/1078) [@falexwolf](https://github.com/falexwolf)
- 👷 Main is the new staging, release the new production branch [PR1077](https://github.com/laminlabs/lamindb/pull/1077) [@falexwolf](https://github.com/falexwolf)
- ♻️ Refactor tutorial [PR1076](https://github.com/laminlabs/lamindb/pull/1076) [@falexwolf](https://github.com/falexwolf)
- ♻️ Add `transform`, `run` & `input_of` to `Dataset` [PR1074](https://github.com/laminlabs/lamindb/pull/1074) [@falexwolf](https://github.com/falexwolf)
- 🚸 Better UX of graphviz installation [PR1075](https://github.com/laminlabs/lamindb/pull/1075) [@sunnyosun](https://github.com/sunnyosun)
- 👷 Run tutorial on anonymous worker [PR1073](https://github.com/laminlabs/lamindb/pull/1073) [@falexwolf](https://github.com/falexwolf)
- ♻️ Use `create_path` [PR1072](https://github.com/laminlabs/lamindb/pull/1072) [@falexwolf](https://github.com/falexwolf)
- 🎨 Make `from_values` synonyms aware, `inspect` bionty aware [PR1066](https://github.com/laminlabs/lamindb/pull/1066) [@sunnyosun](https://github.com/sunnyosun)
- ♻️ Use principled conversion to paths with options inheritance [PR1069](https://github.com/laminlabs/lamindb/pull/1069) [@Koncopd](https://github.com/Koncopd)
- ♻️ Refactor storage [PR1067](https://github.com/laminlabs/lamindb/pull/1067) [@falexwolf](https://github.com/falexwolf)
- 🐛 Fix tree view for anonymous access [PR1065](https://github.com/laminlabs/lamindb/pull/1065) [@falexwolf](https://github.com/falexwolf)
- 🍱 Add the Iris images dev dataset [PR1062](https://github.com/laminlabs/lamindb/pull/1062) [@falexwolf](https://github.com/falexwolf)
- 🚚 `Dataset.from_files` replaced with default constructor [PR1060](https://github.com/laminlabs/lamindb/pull/1060) [@falexwolf](https://github.com/falexwolf)
- 🍱 Added `populate_registries` to some datasets [PR1061](https://github.com/laminlabs/lamindb/pull/1061) [@sunnyosun](https://github.com/sunnyosun)
- ✨ Add `get_feature_set` to `.features` accessor [PR1057](https://github.com/laminlabs/lamindb/pull/1057) [@sunnyosun](https://github.com/sunnyosun)
- 🚸 Allow `reference` & `reference_type` params in `ln.track()` [PR1053](https://github.com/laminlabs/lamindb/pull/1053) [@falexwolf](https://github.com/falexwolf)

## 0.50

### 0.50.1 {small}`2023-08-17`

- ♻️ Refactor `run_context` [PR1030](https://github.com/laminlabs/lamindb/pull/1030) [@falexwolf](https://github.com/falexwolf)
- 🚸 Remove SA as a dependency and use delete dialogue [PR1041](https://github.com/laminlabs/lamindb/pull/1041) [@falexwolf](https://github.com/falexwolf)
- 🚚 Move data lineage to use cases [PR1026](https://github.com/laminlabs/lamindb/pull/1026) [@falexwolf](https://github.com/falexwolf)
- 🚸 Do not allow creating files with none of `key`, `description` and `run` [PR1022](https://github.com/laminlabs/lamindb/pull/1022) [@falexwolf](https://github.com/falexwolf)
- ♻️ Introduce `ParentsAware` [PR1023](https://github.com/laminlabs/lamindb/pull/1023) [@sunnyosun](https://github.com/sunnyosun)
- ♻️ Re-enable file versioning [PR1014](https://github.com/laminlabs/lamindb/pull/1014) [PR1018](https://github.com/laminlabs/lamindb/pull/1018) [@falexwolf](https://github.com/falexwolf)
- 🔊 Fix validate logging [PR1017](https://github.com/laminlabs/lamindb/pull/1017) [@sunnyosun](https://github.com/sunnyosun)
- ✨ Add `InspectResult` [PR1012](https://github.com/laminlabs/lamindb/pull/1012) [@sunnyosun](https://github.com/sunnyosun)
- 🚚 Rename `map_synonyms()` to `standardize()` [PR1046](https://github.com/laminlabs/lamindb/pull/1046) [@sunnyosun](https://github.com/sunnyosun)
- 💄 Improved logging & visualizations [PR1011](https://github.com/laminlabs/lamindb/pull/1011) [@sunnyosun](https://github.com/sunnyosun)

### 0.50.0 {small}`2023-08-07`

- 🚸 Introduce Label validation [PR1006](https://github.com/laminlabs/lamindb/pull/1006) [@falexwolf](https://github.com/falexwolf)
- 🚸 Turn `File.path()` into property [PR1008](https://github.com/laminlabs/lamindb/pull/1008) [@falexwolf](https://github.com/falexwolf)
- 🔥 Remove backward compat for auto-keyed paths [PR1007](https://github.com/laminlabs/lamindb/pull/1007) [@falexwolf](https://github.com/falexwolf)
- ✨ Add `ValidationAware.validate()` [PR1000](https://github.com/laminlabs/lamindb/pull/1000) [@sunnyosun](https://github.com/sunnyosun)
- 🚸 Increase default verbosity to hint-level [PR1005](https://github.com/laminlabs/lamindb/pull/1005) [@falexwolf](https://github.com/falexwolf)
- ♻️ Better treatment of suffixes [PR1004](https://github.com/laminlabs/lamindb/pull/1004) [@falexwolf](https://github.com/falexwolf)
- 🍱 Add nf-core rna-seq test dataset [PR1002](https://github.com/laminlabs/lamindb/pull/1002) [@falexwolf](https://github.com/falexwolf)
- 🚚 Move `add_labels` and `get_labels` from `FeatureManager` to `Data` [PR999](https://github.com/laminlabs/lamindb/pull/999) [@falexwolf](https://github.com/falexwolf)
- 🚚 Rename `FeatureSet.ref_field` to `FeatureSet.registry` [PR997](https://github.com/laminlabs/lamindb/pull/997) [@falexwolf](https://github.com/falexwolf)
- 🚸 Ignore `synonyms` in `from_values` to incentivize explicit curation process [PR993](https://github.com/laminlabs/lamindb/pull/993) [@falexwolf](https://github.com/falexwolf)
- 🚚 Move `tutorial1` one level down [PR989](https://github.com/laminlabs/lamindb/pull/989) [@falexwolf](https://github.com/falexwolf)
- 🚸 Only link validated features [PR987](https://github.com/laminlabs/lamindb/pull/987) [@falexwolf](https://github.com/falexwolf)
- ✅ More fixes and tests for backed [PR988](https://github.com/laminlabs/lamindb/pull/988) [@Koncopd](https://github.com/Koncopd)
- 🚚 Move biology notebooks to `lamin-usecases` [PR985](https://github.com/laminlabs/lamindb/pull/985) [@sunnyosun](https://github.com/sunnyosun)
- 📝 Improve API layout & description [PR986](https://github.com/laminlabs/lamindb/pull/986) [@falexwolf](https://github.com/falexwolf)
- 🚚 Rename `Manager` to `QueryManager` [PR984](https://github.com/laminlabs/lamindb/pull/984) [@falexwolf](https://github.com/falexwolf)
- 🚚 Rename `ORM` to `Registry` [PR982](https://github.com/laminlabs/lamindb/pull/982) [@sunnyosun](https://github.com/sunnyosun)

## 0.49

### 0.49.3 {small}`2023-08-04`

- 🔧 Eliminate index files & upper bound pydantic and other core dependencies by major versions [PR981](https://github.com/laminlabs/lamindb/pull/981) [@falexwolf](https://github.com/falexwolf)

### 0.49.2 {small}`2023-08-03`

- ♻️ Refactor storage root management [PR973](https://github.com/laminlabs/lamindb/pull/973) [@falexwolf](https://github.com/falexwolf)
- 🎨 Make `add_feature_set` acid [PR980](https://github.com/laminlabs/lamindb/pull/980) [@sunnyosun](https://github.com/sunnyosun)
- 📝 Simplify tutorial [PR977](https://github.com/laminlabs/lamindb/pull/977) [@sunnyosun](https://github.com/sunnyosun)
- 🐛 Account for `universal_pathlib` 0.1 release [PR976](https://github.com/laminlabs/lamindb/pull/976) [@Koncopd](https://github.com/Koncopd)

### 0.49.1 {small}`2023-08-01`

- 🎨 Allow passing `QuerySet` to `.features.add_labels()` [PR969](https://github.com/laminlabs/lamindb/pull/969) [@sunnyosun](https://github.com/sunnyosun)
- ⚡️ Speed up & remove auto-generation of feature name [PR963](https://github.com/laminlabs/lamindb/pull/963) [@falexwolf](https://github.com/falexwolf)

### 0.49.0 {small}`2023-07-31`

- ✨ Add a progress monitoring to file.load() and file.stage() [PR959](https://github.com/laminlabs/lamindb/pull/959) [@Koncopd](https://github.com/Koncopd)
- 🚚 Rename `Registry.select()` to `Registry.filter()` [PR958](https://github.com/laminlabs/lamindb/pull/958) [@falexwolf](https://github.com/falexwolf)
- 🚚 Rename `Feature.label_orms` to `Feature.registries` [PR957](https://github.com/laminlabs/lamindb/pull/957) [@falexwolf](https://github.com/falexwolf)
- ♻️ Aggregate ref_field, ref_orm, ref_schema into ref_field [PR955](https://github.com/laminlabs/lamindb/pull/955) [@falexwolf](https://github.com/falexwolf)
- 🎨 Move feature foreignkey from label to link between label and file [PR954](https://github.com/laminlabs/lamindb/pull/954) [@falexwolf](https://github.com/falexwolf)
- 🎨 Enable file/dataset-external feature sets [PR953](https://github.com/laminlabs/lamindb/pull/953) [@falexwolf](https://github.com/falexwolf)

## 0.48 {small}`2023-07-26`

- 📝 Add MuData example [PR945](https://github.com/laminlabs/lamindb/pull/945) [@sunnyosun](https://github.com/sunnyosun)
- 🎨 Add `Modality` [PR940](https://github.com/laminlabs/lamindb/pull/940) [@falexwolf](https://github.com/falexwolf)
- 🎨 Describe based on provenance and features [PR936](https://github.com/laminlabs/lamindb/pull/936) [@sunnyosun](https://github.com/sunnyosun)
- 🚸 Add `File.features` accessor: `FeaturesManager` [PR933](https://github.com/laminlabs/lamindb/pull/933) [@falexwolf](https://github.com/falexwolf)
- 📝 Add spatial data example [PR930](https://github.com/laminlabs/lamindb/pull/930) [@falexwolf](https://github.com/falexwolf)
- ✨ Added `with_children` to `view_parents` [PR921](https://github.com/laminlabs/lamindb/pull/921) [@sunnyosun](https://github.com/sunnyosun)
- 📝 Added comprehensive lineage graph example [PR919](https://github.com/laminlabs/lamindb/pull/919) [@sunnyosun](https://github.com/sunnyosun)
- 🚸 `File` methods work prior to `save()` [PR918](https://github.com/laminlabs/lamindb/pull/918) [@Koncopd](https://github.com/Koncopd)
- 🚚 Integrate `Project` and `Tag` into `Label` [PR914](https://github.com/laminlabs/lamindb/pull/914) [@falexwolf](https://github.com/falexwolf)
- ♻️ Replace `lamin_logger` with `lamin_utils` [PR912](https://github.com/laminlabs/lamindb/pull/912) [@falexwolf](https://github.com/falexwolf)
- ♻️ Replace `Run.name` with `Run.reference_type` and more [PR910](https://github.com/laminlabs/lamindb/pull/910) [@falexwolf](https://github.com/falexwolf)
- ✨ Enable to create `File` outside default storage [PR891](https://github.com/laminlabs/lamindb/pull/891) [@falexwolf](https://github.com/falexwolf)
- ✨ Multi-field search [PR898](https://github.com/laminlabs/lamindb/pull/898) [@sunnyosun](https://github.com/sunnyosun)
- ✨ Better `Feature` management [PR897](https://github.com/laminlabs/lamindb/pull/897) [@falexwolf](https://github.com/falexwolf)
- 🎨 Make `ensembl_gene_id`` unique for `Gene` [PR895](https://github.com/laminlabs/lamindb/pull/895) [@sunnyosun](https://github.com/sunnyosun)
- ⚡ Speed up file creation [PR894](https://github.com/laminlabs/lamindb/pull/894) [@falexwolf](https://github.com/falexwolf)
- 🚚 Rename `Readout` to `ExperimentalFactor` [PR893](https://github.com/laminlabs/lamindb/pull/893) [@sunnyosun](https://github.com/sunnyosun)
- 🎨 `set_abbr` adds `abbr` to `synonyms` [PR892](https://github.com/laminlabs/lamindb/pull/892) [@sunnyosun](https://github.com/sunnyosun)

## 0.47 {small}`2023-07-10`

- ✨ View parents [PR858](https://github.com/laminlabs/lamindb/pull/858) [@sunnyosun](https://github.com/sunnyosun)
- ✨ Track parent notebooks [PR859](https://github.com/laminlabs/lamindb/pull/859) [@falexwolf](https://github.com/falexwolf)
- 🚸 Improve data-lineage guide, auto-track run inputs [PR869](https://github.com/laminlabs/lamindb/pull/869) [@falexwolf](https://github.com/falexwolf)
- 🚸 More reliably detect interactive environments `ln.track()` [PR878](https://github.com/laminlabs/lamindb/pull/878) [@Koncopd](https://github.com/Koncopd)
- 🚸 Case insensitive search [PR877](https://github.com/laminlabs/lamindb/pull/877) [@falexwolf](https://github.com/falexwolf)
- 🚸 Fix bug and add test for tracking multiple parent transforms [PR875](https://github.com/laminlabs/lamindb/pull/875) [@falexwolf](https://github.com/falexwolf)
- ✨ Read remote h5ad in backed mode in `File.from_anndata` [PR871](https://github.com/laminlabs/lamindb/pull/871) [@falexwolf](https://github.com/falexwolf)
- ⚡️ Speed up search [PR868](https://github.com/laminlabs/lamindb/pull/868) [@falexwolf](https://github.com/falexwolf)
- ✨ Added `File.inherit_relationships` [PR867](https://github.com/laminlabs/lamindb/pull/867) [@sunnyosun](https://github.com/sunnyosun)
- ⚡️ Better cache management [PR864](https://github.com/laminlabs/lamindb/pull/864) [@Koncopd](https://github.com/Koncopd)
- ✨ Added `Manager.list()` and allow passing field name [PR863](https://github.com/laminlabs/lamindb/pull/863) [@sunnyosun](https://github.com/sunnyosun)
- ✨ Added `.describe()` for rich repr of related objects [PR862](https://github.com/laminlabs/lamindb/pull/862) [@sunnyosun](https://github.com/sunnyosun)
- ✨ Track `hash_type` [PR861](https://github.com/laminlabs/lamindb/pull/861) [@falexwolf](https://github.com/falexwolf)
- 🚸 It should be possible to use a schema module prior to importing `lamindb` [PR852](https://github.com/laminlabs/lamindb/pull/852) [@falexwolf](https://github.com/falexwolf)
- ✨ Enable `search`, `lookup`, `inspect`, `map_synonyms` from `QuerySet` [PR849](https://github.com/laminlabs/lamindb/pull/849) [@sunnyosun](https://github.com/sunnyosun)

## 0.46 {small}`2023-07-06`

Features:

- ✨ Enable hierarchical metadata, e.g., cell types, tissues, etc. [PR810](https://github.com/laminlabs/lamindb/pull/810) [@sunnyosun](https://github.com/sunnyosun)
- ✨ Add `Dataset` & `Feature` ORMs, e.g., easily track column names of dataframes [PR805](https://github.com/laminlabs/lamindb/pull/805) [@falexwolf](https://github.com/falexwolf)

Breaking changes:

- 🚚 Rename `File.name` to `File.description` [PR824](https://github.com/laminlabs/lamindb/pull/824) [@falexwolf](https://github.com/falexwolf)
- 🚚 Rename `File.featuresets` to `File.feature_sets` [PR805](https://github.com/laminlabs/lamindb/pull/805) [@falexwolf](https://github.com/falexwolf)

UX:

- ✨ Globally set species via `lb.settings.species=` [PR142](https://github.com/laminlabs/lnschema-bionty/pull/142) [falexwolf](https://github.com/falexwolf)
- 🚚 Easy display of many-to-many fields: `QuerySet.df(include=[field__name])` [PR832](https://github.com/laminlabs/lamindb/pull/832) [@falexwolf](https://github.com/falexwolf)
- 🚸 Create new run if notebook is run by different user [PR838](https://github.com/laminlabs/lamindb/pull/838) [@falexwolf](https://github.com/falexwolf)
- 🚸 Speed up bulk saving of records [PR828](https://github.com/laminlabs/lamindb/pull/828) [@sunnyosun](https://github.com/sunnyosun)
- 🚸 Hash large files faster [PR836](https://github.com/laminlabs/lamindb/pull/836) [@falexwolf](https://github.com/falexwolf)
- ✨ Add `from_df` and `from_anndata` to `File` [PR844](https://github.com/laminlabs/lamindb/pull/844) [@falexwolf](https://github.com/falexwolf)
- 🚸 Return locally backed object instead of cloud backed if available [PR840](https://github.com/laminlabs/lamindb/pull/840) [@falexwolf](https://github.com/falexwolf)
- 🚸 Raise more & more user-friendly errors in setup API when instance already setup [PR837](https://github.com/laminlabs/lamindb/pull/837) [@falexwolf](https://github.com/falexwolf)
- 🚸 Better error behavior when no notebook title set in `ln.track()` [PR834](https://github.com/laminlabs/lamindb/pull/834) [@Koncopd](https://github.com/Koncopd)
- 🚸 Store hash for remote files on S3 [PR808](https://github.com/laminlabs/lamindb/pull/808) [@falexwolf](https://github.com/falexwolf)

## 0.45 {small}`2023-06-27`

- ♻️ Replaced `ln.parse` with `Registry.from_values` [PR803](https://github.com/laminlabs/lamindb/pull/803) [@sunnyosun](https://github.com/sunnyosun)
- 🎨 Auto-manage `RunInput` Registry [PR802](https://github.com/laminlabs/lamindb/pull/802) [@falexwolf](https://github.com/falexwolf)

## 0.44

### 0.44.2 {small}`2023-06-23`

- ♻️ Make `zarr` optional [PR800](https://github.com/laminlabs/lamindb/pull/800) [@falexwolf](https://github.com/falexwolf)

### 0.44.1 {small}`2023-06-22`

- ✨ Add `inspect` and `add_synonym` to `Registry` [PR797](https://github.com/laminlabs/lamindb/pull/797) [@sunnyosun](https://github.com/sunnyosun)
- 🔧 Rename extra `nbproject` to `jupyter` and add `fcs` extra to docs [PR798](https://github.com/laminlabs/lamindb/pull/798) [@falexwolf](https://github.com/falexwolf)
- 🚚 Move default storage location from `lndb/` to `.lamindb/` [PR796](https://github.com/laminlabs/lamindb/pull/796) [@falexwolf](https://github.com/falexwolf)
- 🚸 `ln.Folder` becomes `ln.Tag` & directories now modeled as prefixes (as on S3) [PR794](https://github.com/laminlabs/lamindb/pull/794) [@falexwolf](https://github.com/falexwolf)
- ♻️ Refactor storage code [PR792](https://github.com/laminlabs/lamindb/pull/792) [@falexwolf](https://github.com/falexwolf) [@Koncopd](https://github.com/Koncopd)

### 0.44.0 {small}`2023-06-20`

Features:

- 🚸 Idempotency across metadata records & data artifacts [FAQ](https://lamin.ai/docs/faq/idempotency) [PR783](https://github.com/laminlabs/lamindb/pull/783) [@falexwolf](https://github.com/falexwolf)
- ✨ `add_synonym` & `map_synonyms` to enable, e.g., `add_synonym("MyGeneName")` [PR786](https://github.com/laminlabs/lamindb/pull/786) [@sunnyosun](https://github.com/sunnyosun)
- ✨ Backed access for general HDF5 and zarr objects [PR781](https://github.com/laminlabs/lamindb/pull/781) [@Koncopd](https://github.com/Koncopd)

Refactors:

- 🚸 Return records list from `.from_bionty` for multiple matches [PR789](https://github.com/laminlabs/lamindb/pull/789) [@sunnyosun](https://github.com/sunnyosun)
- 🎨 Remove `lnhub-rest` from `lamindb-setup` [PR784](https://github.com/laminlabs/lamindb/pull/784) [bpenteado](https://github.com/bpenteado)
- 🔊 Move logging from stderr to stdout [PR776](https://github.com/laminlabs/lamindb/pull/776) [@falexwolf](https://github.com/falexwolf)

## 0.43.0 {small}`2023-06-15`

Features:

- ✨ Enable `Registry.search()` and improved `Registry.lookup()` [PR771](https://github.com/laminlabs/lamindb/pull/771) [@sunnyosun](https://github.com/sunnyosun)
- 🎨 Consolidate `lnschema_bionty` and upgrade to latest Bionty [PR775](https://github.com/laminlabs/lamindb/pull/775) [@sunnyosun](https://github.com/sunnyosun)
- 🚸 Introduce `ln.settings.storage` to switch default storage [PR773](https://github.com/laminlabs/lamindb/pull/773) [@falexwolf](https://github.com/falexwolf)
- 🚸 Return existing file if hash exists (idempotency) [PR772](https://github.com/laminlabs/lamindb/pull/772) [@falexwolf](https://github.com/falexwolf)
- 🚸 `ln.settings` can now change logging verbosity levels [PR630](https://github.com/laminlabs/lamindb/pull/630) [@falexwolf](https://github.com/falexwolf)

Refactors:

- ♻️ Refactor core schema methods and storage access [PR770](https://github.com/laminlabs/lamindb/pull/770) [@falexwolf](https://github.com/falexwolf)
- 🐛 Make `User.name` nullable again [PR769](https://github.com/laminlabs/lamindb/pull/769) [@falexwolf](https://github.com/falexwolf)
- ✅ Add integrity tests for migrations back [PR768](https://github.com/laminlabs/lamindb/pull/768) [@falexwolf](https://github.com/falexwolf)

## 0.42.0 {small}`2023-06-14`

This is the first release after migrating from SQLModel/SQLAlchemy to Django.

With this, we're hopeful that we get closer to a production-ready 1.0.0 API.

Highlights:

- More robust & simpler automated migrations: `lamin migrate create` & `lamin migrate deploy`
- Simpler query syntax (no joins anymore): `ln.File.select(transform__created_by=user)`
- No need to create a session object to load relationships: access `file.transform` to load a `Transform` object
- No need to write out link models in schemas & generally simplified schema syntax
- Any schema package (`lnschema_myschema`) is now managed as minimal Django app

Breaking changes:

- Renamed `ln.Features` to `ln.FeatureSet` and is now typically instantiated with `FeatureSet.from_iterable()`
- Removed `ln.Session`
- Removed `.join()` (replaced `SelectStmt` with `QuerySet`)
- `.all()` now returns a `QuerySet` and no longer a list (use `.list()` instead)
- Access `Bionty` objects within `lnschema_bionty` via `Registry.bionty()` instead of `Registry.bionty`
- Removed `File.stream()` as all functionality is now provided through `File.backed()`
- Many-to-many fields are now set with `Run.inputs.set()` and appended with `Run.inputs.add()`

Non-breaking changes:

- Vastly simplified dependencies & introduced configurable installation
- Auto-generated storage keys are now of the form `lndb/{id}.{suffix}` rather than just `{id}.{suffix}`
- Renamed `ln.add()` to `ln.save()`
- Introduced `Registry.select()`, `Registry.save()`, and `Registry.delete()`
- Better tracking & linking of Bionty sources in `lnschema_bionty`

Additional notes:

- Consolidated docs and auto-generate upon push events to lamindb main
- Consolidated submodules (renamed `lndb` to `lamindb-setup`, removed `lndb-storage`)

The main downsides of migrating to Django are:

- Currently only one LaminDB instance loadable per Python session
- Type hints & constructor signatures are less pythonic (SQLModel uses less magic than Django) and lead to idiosyncrasies in model definition (nullable defaults) and validation (validation at the Registry-level is more manual as Django foresees validation at the Form level)
- SQLAlchemy provides the more powerful Registry, and there might be future use cases that will require them

Complete list of changes below.

<!-- prettier-ignore -->
Name | PR | Developer | Date | Version
--- | --- | --- | --- | ---
♻️ Use `TransformType` | [763](https://github.com/laminlabs/lamindb/pull/763) | [falexwolf](https://github.com/falexwolf) | 2023-06-13 |
👷 Dispatch to lamin-examples & redun-lamin-fasta | [762](https://github.com/laminlabs/lamindb/pull/762) | [falexwolf](https://github.com/falexwolf) | 2023-06-13 |
🔥 Remove `File.stream()` | [761](https://github.com/laminlabs/lamindb/pull/761) | [falexwolf](https://github.com/falexwolf) | 2023-06-13 | 0.42a9
🚸 Prefix auto-storage-key with `lndb/` | [757](https://github.com/laminlabs/lamindb/pull/757) | [falexwolf](https://github.com/falexwolf) | 2023-06-12 |
✨ Delete storage in `File.delete()` | [754](https://github.com/laminlabs/lamindb/pull/754) | [Koncopd](https://github.com/Koncopd) | 2023-06-12 |
✅ Add more tests for File init | [755](https://github.com/laminlabs/lamindb/pull/755) | [falexwolf](https://github.com/falexwolf) | 2023-06-12 |
💚 Remove test paths from pyproject.toml | [753](https://github.com/laminlabs/lamindb/pull/753) | [falexwolf](https://github.com/falexwolf) | 2023-06-11 |
✨ Add `to_adata()` method to `AnnDataAccessor` | [752](https://github.com/laminlabs/lamindb/pull/752) | [Koncopd](https://github.com/Koncopd) | 2023-06-11 |
🚚 Rename `Featureset` to `FeatureSet` | [750](https://github.com/laminlabs/lamindb/pull/750) | [sunnyosun](https://github.com/sunnyosun) | 2023-06-10 |
⚠️ Refactor save - it no longer returns records | [742](https://github.com/laminlabs/lamindb/pull/742) | [falexwolf](https://github.com/falexwolf) | 2023-06-10 |
📝 Re-organize biology guides | [740](https://github.com/laminlabs/lamindb/pull/740) | [falexwolf](https://github.com/falexwolf) | 2023-06-10 |
✨ Populate `bionty_version` in `ln.parse` | [739](https://github.com/laminlabs/lamindb/pull/739) | [sunnyosun](https://github.com/sunnyosun) | 2023-06-09 | 0.42a8
⚡ Improved multifield query in `ln.parse` | [736](https://github.com/laminlabs/lamindb/pull/736) | [sunnyosun](https://github.com/sunnyosun) | 2023-06-09 |
⬆️ First stable Django release of `lnschema-bionty` | [733](https://github.com/laminlabs/lamindb/pull/733) | [sunnyosun](https://github.com/sunnyosun) | 2023-06-09 |
🚸 Validate required fields | [735](https://github.com/laminlabs/lamindb/pull/735) | [falexwolf](https://github.com/falexwolf) | 2023-06-09 |
📝 Integrate `lnschema_bionty` into reference | [732](https://github.com/laminlabs/lamindb/pull/732) | [falexwolf](https://github.com/falexwolf) | 2023-06-08 | 0.42a7
🚸 Add `select` method to `BaseORM` | [730](https://github.com/laminlabs/lamindb/pull/730) | [falexwolf](https://github.com/falexwolf) | 2023-06-08 |
📝 Overhaul README | [728](https://github.com/laminlabs/lamindb/pull/728) | [falexwolf](https://github.com/falexwolf) | 2023-06-08 |
♻️ Refactored features code | [731](https://github.com/laminlabs/lamindb/pull/731) | [sunnyosun](https://github.com/sunnyosun) | 2023-06-08 |
📝 Add `configure` guide instead of setup | [727](https://github.com/laminlabs/lamindb/pull/727) | [falexwolf](https://github.com/falexwolf) | 2023-06-08 |
🧪 Add tests for `folder.tree()` | [726](https://github.com/laminlabs/lamindb/pull/726) | [sunnyosun](https://github.com/sunnyosun) | 2023-06-08 |
🚚 Renamed `BiontyVersions` to `BiontySource` | [725](https://github.com/laminlabs/lamindb/pull/725) | [sunnyosun](https://github.com/sunnyosun) | 2023-06-07 |
🔊 `folder.tree` can only be used with existing folders in storage | [724](https://github.com/laminlabs/lamindb/pull/724) | [sunnyosun](https://github.com/sunnyosun) | 2023-06-07 |
👷 Coverage in separate job | [722](https://github.com/laminlabs/lamindb/pull/722) | [falexwolf](https://github.com/falexwolf) | 2023-06-07 |
🎨 Import ORMs from .models before reload | [723](https://github.com/laminlabs/lamindb/pull/723) | [sunnyosun](https://github.com/sunnyosun) | 2023-06-07 |
🔥 Remove `ln.nb` | [721](https://github.com/laminlabs/lamindb/pull/721) | [falexwolf](https://github.com/falexwolf) | 2023-06-07 |
🏗️ Re-architect transform id | [720](https://github.com/laminlabs/lamindb/pull/720) | [falexwolf](https://github.com/falexwolf) | 2023-06-06 | 0.42a6
🚚 Rename `ln.add` to `ln.save` | [719](https://github.com/laminlabs/lamindb/pull/719) | [falexwolf](https://github.com/falexwolf) | 2023-06-05 | 0.42a5
⬆️ Upgrade lnschema-core to 0.35a5 | [718](https://github.com/laminlabs/lamindb/pull/718) | [sunnyosun](https://github.com/sunnyosun) | 2023-06-05 |
🚚 Migrate `lnschema-bionty` to Django | [716](https://github.com/laminlabs/lamindb/pull/716) | [falexwolf](https://github.com/falexwolf) | 2023-06-05 | 0.42a2
♻️ Polish core schema | [717](https://github.com/laminlabs/lamindb/pull/717) | [falexwolf](https://github.com/falexwolf) | 2023-06-05 |
🐛 Fix delete for File | [715](https://github.com/laminlabs/lamindb/pull/715) | [Koncopd](https://github.com/Koncopd) | 2023-06-04 |
📝 Rework the stream notebook | [714](https://github.com/laminlabs/lamindb/pull/714) | [Koncopd](https://github.com/Koncopd) | 2023-06-04 |
💚 Point to lamindb-setup main & fix session | [713](https://github.com/laminlabs/lamindb/pull/713) | [falexwolf](https://github.com/falexwolf) | 2023-06-04 |
➖ Move nbproject to extra dependencies | [711](https://github.com/laminlabs/lamindb/pull/711) | [Koncopd](https://github.com/Koncopd) | 2023-06-04 |
🔥 Delete SQLAlchemy related content | [710](https://github.com/laminlabs/lamindb/pull/710) | [falexwolf](https://github.com/falexwolf) | 2023-06-04 |
➕ Pin boto3 in aws | [712](https://github.com/laminlabs/lamindb/pull/712) | [Koncopd](https://github.com/Koncopd) | 2023-06-04 |
🔥 Remove SQLAlchemy tests | [709](https://github.com/laminlabs/lamindb/pull/709) | [falexwolf](https://github.com/falexwolf) | 2023-06-04 |
♻️ Absorb `DjangoORM.create()` in `DjangoORM.__init__()` | [707](https://github.com/laminlabs/lamindb/pull/707) | [falexwolf](https://github.com/falexwolf) | 2023-06-03 |
🐛 Disentangle keys in storage related test notebooks | [708](https://github.com/laminlabs/lamindb/pull/708) | [Koncopd](https://github.com/Koncopd) | 2023-06-03 |
🚸 ln.track improvements | [704](https://github.com/laminlabs/lamindb/pull/704) | [Koncopd](https://github.com/Koncopd) | 2023-06-03 |
🏗️ Enable Django backend (part 2) | [702](https://github.com/laminlabs/lamindb/pull/702) | [falexwolf](https://github.com/falexwolf) | 2023-06-02 |
🎨 Simplified track sample-level metadata | [705](https://github.com/laminlabs/lamindb/pull/705) | [sunnyosun](https://github.com/sunnyosun) | 2023-06-02 |
🔊 Add more loggings to `ln.parse` | [703](https://github.com/laminlabs/lamindb/pull/703) | [sunnyosun](https://github.com/sunnyosun) | 2023-06-02 | 0.42a1
♻️ Refactored feature parsing and `ln.parse` | [701](https://github.com/laminlabs/lamindb/pull/701) | [sunnyosun](https://github.com/sunnyosun) | 2023-06-01 |
🚚 Move lndb-storage back into lamindb | [700](https://github.com/laminlabs/lamindb/pull/700) | [falexwolf](https://github.com/falexwolf) | 2023-06-01 |
🚚 Rename `lndb` to `lamindb_setup` | [699](https://github.com/laminlabs/lamindb/pull/699) | [falexwolf](https://github.com/falexwolf) | 2023-06-01 |
🏗️ Add Django backend (setup) | [697](https://github.com/laminlabs/lamindb/pull/697) | [falexwolf](https://github.com/falexwolf) | 2023-05-31 |
⬆️ Update lndb to 0.45.0 | [698](https://github.com/laminlabs/lamindb/pull/698) | [Koncopd](https://github.com/Koncopd) | 2023-05-31 | 0.41.2
⬆️ Upgrade lnschema-bionty | [696](https://github.com/laminlabs/lamindb/pull/696) | [sunnyosun](https://github.com/sunnyosun) | 2023-05-30 | 0.41.1
🚑 Fix species config | [695](https://github.com/laminlabs/lamindb/pull/695) | [sunnyosun](https://github.com/sunnyosun) | 2023-05-28 | 0.41.0
🎨 Clean up CI more | [694](https://github.com/laminlabs/lamindb/pull/694) | [falexwolf](https://github.com/falexwolf) | 2023-05-28 |
➖ Remove scanpy as test dependency | [693](https://github.com/laminlabs/lamindb/pull/693) | [falexwolf](https://github.com/falexwolf) | 2023-05-28 |
👷 Fix coverage for lndb-storage | [692](https://github.com/laminlabs/lamindb/pull/692) | [falexwolf](https://github.com/falexwolf) | 2023-05-28 |
➖ Do not install storage extras by default | [691](https://github.com/laminlabs/lamindb/pull/691) | [Koncopd](https://github.com/Koncopd) | 2023-05-28 | 0.41a4
👷 Bring back nox session | [690](https://github.com/laminlabs/lamindb/pull/690) | [falexwolf](https://github.com/falexwolf) | 2023-05-28 |
📝 Refactor guide notebooks | [689](https://github.com/laminlabs/lamindb/pull/689) | [falexwolf](https://github.com/falexwolf) | 2023-05-28 |
✨ Add `DataFrame` support for `ln.parse` | [688](https://github.com/laminlabs/lamindb/pull/688) | [sunnyosun](https://github.com/sunnyosun) | 2023-05-28 |
👷 Refactor tests | [687](https://github.com/laminlabs/lamindb/pull/687) | [falexwolf](https://github.com/falexwolf) | 2023-05-27 |
✨ Subsettable backed `AnnData` | [668](https://github.com/laminlabs/lamindb/pull/668) | [Koncopd](https://github.com/Koncopd) | 2023-05-27 |
📝 Remove setup notebook | [686](https://github.com/laminlabs/lamindb/pull/686) | [falexwolf](https://github.com/falexwolf) | 2023-05-27 |
🏗️ Remove SQL-level schema modules | [685](https://github.com/laminlabs/lamindb/pull/685) | [falexwolf](https://github.com/falexwolf) | 2023-05-26 | 0.41a3
⬆️ Upgrade lnschema-bionty to 0.17.1 | [684](https://github.com/laminlabs/lamindb/pull/684) | [sunnyosun](https://github.com/sunnyosun) | 2023-05-25 | 0.40.7
⬇️ Downgrade lnschema-bionty to 0.16.5 | [683](https://github.com/laminlabs/lamindb/pull/683) | [sunnyosun](https://github.com/sunnyosun) | 2023-05-25 | 0.40.6
♻️ Refactor types | [681](https://github.com/laminlabs/lamindb/pull/681) | [falexwolf](https://github.com/falexwolf) | 2023-05-23 | 0.40.5
♻️ Refactor `BaseORM` | [679](https://github.com/laminlabs/lamindb/pull/679) | [falexwolf](https://github.com/falexwolf) | 2023-05-17 |
🚸 Pre-join some cheap relationships | [678](https://github.com/laminlabs/lamindb/pull/678) | [falexwolf](https://github.com/falexwolf) | 2023-05-16 |
📝 Improve wording | [677](https://github.com/laminlabs/lamindb/pull/677) | [Zethson](https://github.com/Zethson) | 2023-05-16 |
✨ Added Treatment table | [675](https://github.com/laminlabs/lamindb/pull/675) | [sunnyosun](https://github.com/sunnyosun) | 2023-05-15 | 0.40.3
🚑 Fix ln.Features table name for postgres | [674](https://github.com/laminlabs/lamindb/pull/674) | [sunnyosun](https://github.com/sunnyosun) | 2023-05-11 | 0.40.1
📝 Update ontology guide | [673](https://github.com/laminlabs/lamindb/pull/673) | [falexwolf](https://github.com/falexwolf) | 2023-05-11 |
🎨 Deprecated `data` in `ln.Features`, replace with `iterable` | [672](https://github.com/laminlabs/lamindb/pull/672) | [sunnyosun](https://github.com/sunnyosun) | 2023-05-11 |
🎨 Replace reference with field for `ln.Features` <span class="badge badge-warning">Breaking</span> | [671](https://github.com/laminlabs/lamindb/pull/671) | [sunnyosun](https://github.com/sunnyosun) | 2023-05-09 | 0.40.0
⬆️ Update bionty to 0.13 | [670](https://github.com/laminlabs/lamindb/pull/670) | [sunnyosun](https://github.com/sunnyosun) | 2023-05-09 | 0.39.8
📝 Polish | [667](https://github.com/laminlabs/lamindb/pull/667) | [falexwolf](https://github.com/falexwolf) | 2023-04-28 |
📝 Polish docs | [666](https://github.com/laminlabs/lamindb/pull/666) | [falexwolf](https://github.com/falexwolf) | 2023-04-28 |
⬆️ Update to lndb 0.44.7 | [665](https://github.com/laminlabs/lamindb/pull/665) | [sunnyosun](https://github.com/sunnyosun) | 2023-04-28 | 0.39.7
🚸 Do not require session for `is_run_input` | [664](https://github.com/laminlabs/lamindb/pull/664) | [falexwolf](https://github.com/falexwolf) | 2023-04-28 | 0.39.6
⬆️ Allow load with storage | [663](https://github.com/laminlabs/lamindb/pull/663) | [sunnyosun](https://github.com/sunnyosun) | 2023-04-27 | 0.39.5
💥 Switch to methods bionty.df(), bionty.lookup() <span class="badge badge-warning">Breaking</span> | [662](https://github.com/laminlabs/lamindb/pull/662) | [sunnyosun](https://github.com/sunnyosun) | 2023-04-27 | 0.39.4
⬆️ Upgrade lndb | [661](https://github.com/laminlabs/lamindb/pull/661) | [fredericenard](https://github.com/fredericenard) | 2023-04-26 | 0.39.3
✨ Enable database entries lookup | [660](https://github.com/laminlabs/lamindb/pull/660) | [sunnyosun](https://github.com/sunnyosun) | 2023-04-24 | 0.39.2
✨ Introduce `File.backed()` | [659](https://github.com/laminlabs/lamindb/pull/659) | [falexwolf](https://github.com/falexwolf) | 2023-04-24 |
✨ Introduced `ln.parse()` | [658](https://github.com/laminlabs/lamindb/pull/658) | [sunnyosun](https://github.com/sunnyosun) | 2023-04-24 | 0.39.1
🎨 Refactor `ln.File` | [657](https://github.com/laminlabs/lamindb/pull/657) | [falexwolf](https://github.com/falexwolf) | 2023-04-24 | 0.39.0
✨ Inherit fsspec kwargs from root and move in root check | [655](https://github.com/laminlabs/lamindb/pull/655) | [Koncopd](https://github.com/Koncopd) | 2023-04-23 |
📌 Check anndata version in `File.subset` | [commit](https://github.com/laminlabs/lamindb/commit/06652591321665f0bed970f36e5e449d3d430945) | [Koncopd](https://github.com/Koncopd) | 2023-04-23 |
🐛 Fix VS Code notebook reinitialization in track | [654](https://github.com/laminlabs/lamindb/pull/654) | [Koncopd](https://github.com/Koncopd) | 2023-04-23 |
🚸 Add `File.subset` | [653](https://github.com/laminlabs/lamindb/pull/653) | [falexwolf](https://github.com/falexwolf) | 2023-04-22 | 0.39rc1
🚚 Replace `lnschema-wetlab` with `lnbase-biolab` and `lnschema-lamin1` | [651](https://github.com/laminlabs/lamindb/pull/651) | [falexwolf](https://github.com/falexwolf) | 2023-04-22 |
📝 Add an export example for `ln.schema.view()` | [649](https://github.com/laminlabs/lamindb/pull/649) | [falexwolf](https://github.com/falexwolf) | 2023-04-21 |
🐛 Fix population of `transform_id` in `File` in edge cases | [648](https://github.com/laminlabs/lamindb/pull/648) | [falexwolf](https://github.com/falexwolf) | 2023-04-21 | 0.38.3
🚸 Allow registering local postgres instances on the hub | [647](https://github.com/laminlabs/lamindb/pull/647) | [falexwolf](https://github.com/falexwolf) | 2023-04-21 | 0.38.2
🚸 Improve error message for notebook tracking | [643](https://github.com/laminlabs/lamindb/pull/643) | [falexwolf](https://github.com/falexwolf) | 2023-04-19 |
van -> can | [641](https://github.com/laminlabs/lamindb/pull/641) | [ThomVett](https://github.com/ThomVett) | 2023-04-19 |
⚡ Improved feature parsing speed | [640](https://github.com/laminlabs/lamindb/pull/640) | [sunnyosun](https://github.com/sunnyosun) | 2023-04-19 | 0.38.1
⬆️ Compatibility with new hub | [639](https://github.com/laminlabs/lamindb/pull/639) | [falexwolf](https://github.com/falexwolf) | 2023-04-18 | 0.38.0
💥 New calling patterns for lnschema-bionty | [633](https://github.com/laminlabs/lamindb/pull/633) | [falexwolf](https://github.com/falexwolf) | 2023-04-18 |
✅ Use nbproject-test directly | [638](https://github.com/laminlabs/lamindb/pull/638) | [Koncopd](https://github.com/Koncopd) | 2023-04-18 |
🚸 Use relative path in key | [636](https://github.com/laminlabs/lamindb/pull/636) | [falexwolf](https://github.com/falexwolf) | 2023-04-18 | 0.37.2
🩹 Unpack `notebook_path` correctly | [632](https://github.com/laminlabs/lamindb/pull/632) | [Koncopd](https://github.com/Koncopd) | 2023-04-17 | 0.37.1
🚚 Better names, more relationships directly on `File` <span class="badge badge-warning">Breaking</span> | [631](https://github.com/laminlabs/lamindb/pull/631) | [falexwolf](https://github.com/falexwolf) | 2023-04-16 | 0.37.0
✨ Add Google Colab integration | [628](https://github.com/laminlabs/lamindb/pull/628) | [falexwolf](https://github.com/falexwolf) | 2023-04-16 |
🚸 Improve notebook tracking UX | [627](https://github.com/laminlabs/lamindb/pull/627) | [falexwolf](https://github.com/falexwolf) | 2023-04-14 |
🎨 Simplify `ln.track()` and add `app` transform type | [624](https://github.com/laminlabs/lamindb/pull/624) | [falexwolf](https://github.com/falexwolf) | 2023-04-12 | 0.36.3
🐛 Fix initialization of new notebooks | [623](https://github.com/laminlabs/lamindb/pull/623) | [Koncopd](https://github.com/Koncopd) | 2023-04-10 | 0.36.1
🚸 Filename in `File.name`, new `File.key` and `Folder.key`, robustness overhaul <span class="badge badge-warning">Breaking</span> | [614](https://github.com/laminlabs/lamindb/pull/614) | [falexwolf](https://github.com/falexwolf) | 2023-04-08 | 0.36.0
✨ Introduce `File.stage()` and `File.replace()` | [611](https://github.com/laminlabs/lamindb/pull/611) | [Koncopd](https://github.com/Koncopd) | 2023-04-03 | 0.35.6
🚸 More robust ontology version tracking | [605](https://github.com/laminlabs/lamindb/pull/605) | [Zethson](https://github.com/Zethson) | 2023-04-01 |
🚚 Import `lnschema_bionty` instead of `bionty` | [591](https://github.com/laminlabs/lamindb/pull/591) | [sunnyosun](https://github.com/sunnyosun) | 2023-03-28 | 0.35.5
🚚 Move all core entities to root level | [590](https://github.com/laminlabs/lamindb/pull/590) | [falexwolf](https://github.com/falexwolf) | 2023-03-27 | 0.35.4
🚸 Polish guide | [589](https://github.com/laminlabs/lamindb/pull/589) | [falexwolf](https://github.com/falexwolf) | 2023-03-27 | 0.35.3
🚚 Move `Readout` from `lnschema-wetlab` to `lnschema-bionty` | [588](https://github.com/laminlabs/lamindb/pull/588) | [sunnyosun](https://github.com/sunnyosun) | 2023-03-27 | 0.35.2
♻️ Fix `parsing_id` | [587](https://github.com/laminlabs/lamindb/pull/587) | [sunnyosun](https://github.com/sunnyosun) | 2023-03-27 | 0.35.1
🚚 Rename `DObject` to `File` and `DFolder` to `Folder` <span class="badge badge-warning">Breaking</span> | [586](https://github.com/laminlabs/lamindb/pull/586) | [falexwolf](https://github.com/falexwolf) | 2023-03-25 | 0.35.0
🚚 Introduce `ln.track()` to replace `ln.nb.header()` | [585](https://github.com/laminlabs/lamindb/pull/585) | [falexwolf](https://github.com/falexwolf) | 2023-03-24 |
🏗️ Combine `Notebook` and `Pipeline` into `Transform` <span class="badge badge-warning">Breaking</span> | [584](https://github.com/laminlabs/lamindb/pull/584) | [falexwolf](https://github.com/falexwolf) | 2023-03-23 |
🚑 Fix optional dependencies | [583](https://github.com/laminlabs/lamindb/pull/583) | [sunnyosun](https://github.com/sunnyosun) | 2023-03-22 | 0.34.2
🔥 Remove `Usage` | [581](https://github.com/laminlabs/lamindb/pull/581) | [falexwolf](https://github.com/falexwolf) | 2023-03-22 |
⬆️ Updated `CellMarker` asset | [578](https://github.com/laminlabs/lamindb/pull/578) | [sunnyosun](https://github.com/sunnyosun) | 2023-03-21 | 0.34.1
⬆️ Upgrade to `bionty` 0.9 | [575](https://github.com/laminlabs/lamindb/pull/575) | [sunnyosun](https://github.com/sunnyosun) | 2023-03-20 | 0.34.0
♻️ Move storage-related code to `lndb-storage` | [560](https://github.com/laminlabs/lamindb/pull/560) | [Koncopd](https://github.com/Koncopd) | 2023-03-20 |
🚑 Fix gene id | [573](https://github.com/laminlabs/lamindb/pull/573) | [sunnyosun](https://github.com/sunnyosun) | 2023-03-17 |
🐛 Check if env variable exists before trying to get his value | [572](https://github.com/laminlabs/lamindb/pull/572) | [fredericenard](https://github.com/fredericenard) | 2023-03-15 |
🚸 Do not yet show deprecation warning for `ln.nb.header()` | [commit](https://github.com/laminlabs/lamindb/commit/2bb80c546434e37b2fa2b0c5b38c92a379aff793) | [falexwolf](https://github.com/falexwolf) | 2023-03-15 | 0.33.4
👷 Restore streaming test from cloud | [569](https://github.com/laminlabs/lamindb/pull/569) | [sunnyosun](https://github.com/sunnyosun) | 2023-03-15 | 0.33.3
⬇️ typeguard<3.0.0 | [568](https://github.com/laminlabs/lamindb/pull/568) | [sunnyosun](https://github.com/sunnyosun) | 2023-03-15 |
➖ Remove typeguard dependency | [567](https://github.com/laminlabs/lamindb/pull/567) | [sunnyosun](https://github.com/sunnyosun) | 2023-03-15 | 0.33.2
⬆️ Update core | [566](https://github.com/laminlabs/lamindb/pull/566) | [sunnyosun](https://github.com/sunnyosun) | 2023-03-15 | 0.33.1
🚚 Replace `ln.nb.header()` with `ln.Run()` except in `faq/nb` | [564](https://github.com/laminlabs/lamindb/pull/564) | [falexwolf](https://github.com/falexwolf) | 2023-03-14 | 0.33.0
🚸 Smart about `global_context` and `load_latest` when run from notebook | [563](https://github.com/laminlabs/lamindb/pull/563) | [falexwolf](https://github.com/falexwolf) | 2023-03-14 |
✨ `ln.Features` | [562](https://github.com/laminlabs/lamindb/pull/562) | [sunnyosun](https://github.com/sunnyosun) | 2023-03-14 |
🏗️ Introduce `lamindb.context` and enable `ln.Run` to create contexts | [561](https://github.com/laminlabs/lamindb/pull/561) | [falexwolf](https://github.com/falexwolf) | 2023-03-13 |
📝 Improve the docstrings of `ln.save` and `ln.delete` | [559](https://github.com/laminlabs/lamindb/pull/559) | [sunnyosun](https://github.com/sunnyosun) | 2023-03-10 |
📝 Hide CI related cells in notebooks | [558](https://github.com/laminlabs/lamindb/pull/558) | [sunnyosun](https://github.com/sunnyosun) | 2023-03-10 |
📝 Update docs to clarify sign up and log in | [557](https://github.com/laminlabs/lamindb/pull/557) | [lawrlee](https://github.com/lawrlee) | 2023-03-10 |
📝 Prettier species query | [555](https://github.com/laminlabs/lamindb/pull/555) | [falexwolf](https://github.com/falexwolf) | 2023-03-09 | 0.32.0
📝 Refactor docs sidebar | [553](https://github.com/laminlabs/lamindb/pull/553) | [sunnyosun](https://github.com/sunnyosun) | 2023-03-09 |
⬆️ Upgrade `ln.setup` | [554](https://github.com/laminlabs/lamindb/pull/554) | [falexwolf](https://github.com/falexwolf) | 2023-03-09 |
🔥 Remove `ln.knowledge` <span class="badge badge-warning">Breaking</span> | [552](https://github.com/laminlabs/lamindb/pull/552) | [sunnyosun](https://github.com/sunnyosun) | 2023-03-09 |
➖ Remove bionty as a dependency | [551](https://github.com/laminlabs/lamindb/pull/551) | [sunnyosun](https://github.com/sunnyosun) | 2023-03-09 | 0.32.0rc1
📝 Replace `ln.knowledge` with `bionty` in docs | [547](https://github.com/laminlabs/lamindb/pull/547) | [falexwolf](https://github.com/falexwolf) | 2023-03-07 | 0.31.1
📝 Link FAQ and guide to session and notebook API | [550](https://github.com/laminlabs/lamindb/pull/550) | [falexwolf](https://github.com/falexwolf) | 2023-03-07 |
⬆️ Upgrade lndb | [549](https://github.com/laminlabs/lamindb/pull/549) | [fredericenard](https://github.com/fredericenard) | 2023-03-07 | 0.31.0
⬆️ Upgrade lndb | [548](https://github.com/laminlabs/lamindb/pull/548) | [fredericenard](https://github.com/fredericenard) | 2023-03-07 |
📝 Move ingest-folder back to faq | [545](https://github.com/laminlabs/lamindb/pull/545) | [falexwolf](https://github.com/falexwolf) | 2023-03-06 |
🚸 Improve clarity of `no_source` error message | [543](https://github.com/laminlabs/lamindb/pull/543) | [falexwolf](https://github.com/falexwolf) | 2023-03-06 |
💚 Fix CI | [542](https://github.com/laminlabs/lamindb/pull/542) | [sunnyosun](https://github.com/sunnyosun) | 2023-03-06 |
🚚 Rename data objects to data or datasets in titles | [541](https://github.com/laminlabs/lamindb/pull/541) | [sunnyosun](https://github.com/sunnyosun) | 2023-03-06 |
📝 Duplicate README to the guide landing page | [540](https://github.com/laminlabs/lamindb/pull/540) | [sunnyosun](https://github.com/sunnyosun) | 2023-03-06 |
✨ Allow uploading zarr from local paths | [539](https://github.com/laminlabs/lamindb/pull/539) | [Koncopd](https://github.com/Koncopd) | 2023-03-05 |
📝 Prettify guide landing page | [537](https://github.com/laminlabs/lamindb/pull/537) | [sunnyosun](https://github.com/sunnyosun) | 2023-03-05 |
🚑 Fix upsert for dobject and CI | [535](https://github.com/laminlabs/lamindb/pull/535) | [sunnyosun](https://github.com/sunnyosun) | 2023-03-03 | 0.30.3
💄 Simplify docs | [534](https://github.com/laminlabs/lamindb/pull/534) | [falexwolf](https://github.com/falexwolf) | 2023-03-03 |
🚸 Do not error, just warn upon installation of lamin | [533](https://github.com/laminlabs/lamindb/pull/533) | [falexwolf](https://github.com/falexwolf) | 2023-03-02 | 0.30.2
🚑 Fix iterdir to not list itself | [532](https://github.com/laminlabs/lamindb/pull/532) | [sunnyosun](https://github.com/sunnyosun) | 2023-03-02 |
🔥 Remove `lamin` dependency | [531](https://github.com/laminlabs/lamindb/pull/531) | [falexwolf](https://github.com/falexwolf) | 2023-03-01 | 0.30.1
🚑 Fix for listing cloud dir | [528](https://github.com/laminlabs/lamindb/pull/528) | [sunnyosun](https://github.com/sunnyosun) | 2023-03-01 | 0.30.0
🎨 Add `lamindb.setup` API as an alternative for the CLI | [530](https://github.com/laminlabs/lamindb/pull/530) | [falexwolf](https://github.com/falexwolf) | 2023-03-01 |
🚚 Rename `dfolder.get_dobject` to `get`, allow passing a subdirectory | [527](https://github.com/laminlabs/lamindb/pull/527) | [sunnyosun](https://github.com/sunnyosun) | 2023-03-01 |
🚑 Ensure filepath is absolute | [526](https://github.com/laminlabs/lamindb/pull/526) | [sunnyosun](https://github.com/sunnyosun) | 2023-03-01 |
🎨 Allow passing a list of relpaths to `dfolder.get_dobject` | [524](https://github.com/laminlabs/lamindb/pull/524) | [sunnyosun](https://github.com/sunnyosun) | 2023-02-28 |
✨ Infer filesystem for anndata read and write | [522](https://github.com/laminlabs/lamindb/pull/522) | [Koncopd](https://github.com/Koncopd) | 2023-02-25 |
📝 Replace `lndb` with `lamin` in docs | [521](https://github.com/laminlabs/lamindb/pull/521) | [falexwolf](https://github.com/falexwolf) | 2023-02-25 | 0.29.1
📝 Add nbproject note box to ingest guide | [520](https://github.com/laminlabs/lamindb/pull/520) | [bpenteado](https://github.com/bpenteado) | 2023-02-24 |
🍱 `DFolder.get_dobject` | [519](https://github.com/laminlabs/lamindb/pull/519) | [sunnyosun](https://github.com/sunnyosun) | 2023-02-24 |
⬆️  Upgrade lnschema-core to 0.28.6 | [517](https://github.com/laminlabs/lamindb/pull/517) | [bpenteado](https://github.com/bpenteado) | 2023-02-23 | 0.29.0
💥 Move `lns.DFolder` to `ln.DFolder` | [510](https://github.com/laminlabs/lamindb/pull/510) | [sunnyosun](https://github.com/sunnyosun) | 2023-02-23 |
🐛 Fix trailing slash within lamindb | [516](https://github.com/laminlabs/lamindb/pull/516) | [falexwolf](https://github.com/falexwolf) | 2023-02-22 | 0.28.5
👷 Fix CI config | [515](https://github.com/laminlabs/lamindb/pull/515) | [falexwolf](https://github.com/falexwolf) | 2023-02-22 | 0.28.4
🐛 Fix version check | [514](https://github.com/laminlabs/lamindb/pull/514) | [falexwolf](https://github.com/falexwolf) | 2023-02-22 | 0.28.3
👷 Collect docs artifacts | [513](https://github.com/laminlabs/lamindb/pull/513) | [falexwolf](https://github.com/falexwolf) | 2023-02-22 | 0.28.2
🐛 Fix trailing slash in storage root | [512](https://github.com/laminlabs/lamindb/pull/512) | [falexwolf](https://github.com/falexwolf) | 2023-02-22 | 0.28.1
🚚 Rename `DObject.run_id` to `DObject.source_id` | [509](https://github.com/laminlabs/lamindb/pull/509) | [falexwolf](https://github.com/falexwolf) | 2023-02-21 | 0.28.0
🐛 Another occurance of local filepath | [508](https://github.com/laminlabs/lamindb/pull/508) | [falexwolf](https://github.com/falexwolf) | 2023-02-21 |
👷 Better CI and better arg validation | [505](https://github.com/laminlabs/lamindb/pull/505) | [falexwolf](https://github.com/falexwolf) | 2023-02-21 |
🐛 Fix tracking local existing data | [506](https://github.com/laminlabs/lamindb/pull/506) | [falexwolf](https://github.com/falexwolf) | 2023-02-21 |
🚑 Fix parents for existing data ingestion | [504](https://github.com/laminlabs/lamindb/pull/504) | [sunnyosun](https://github.com/sunnyosun) | 2023-02-21 |
🎨 Disable multiple select results for .add and .delete by fields | [502](https://github.com/laminlabs/lamindb/pull/502) | [sunnyosun](https://github.com/sunnyosun) | 2023-02-20 | 0.28rc1
⚡️ Replace CloudPath with UPath | [501](https://github.com/laminlabs/lamindb/pull/501) | [Koncopd](https://github.com/Koncopd) | 2023-02-19 |
👷 Move fixtures from nox to conftest | [500](https://github.com/laminlabs/lamindb/pull/500) | [falexwolf](https://github.com/falexwolf) | 2023-02-17 |
📝 Simplify output for dfolder faq | [499](https://github.com/laminlabs/lamindb/pull/499) | [sunnyosun](https://github.com/sunnyosun) | 2023-02-17 |
📝 Guide to ingest a folder | [496](https://github.com/laminlabs/lamindb/pull/496) | [sunnyosun](https://github.com/sunnyosun) | 2023-02-17 |
🚑 Fix determination of sqlite vs postgres | [497](https://github.com/laminlabs/lamindb/pull/497) | [falexwolf](https://github.com/falexwolf) | 2023-02-17 | 0.27.2
🎨 Improve `lndb` (lamindb manager) architecture | [495](https://github.com/laminlabs/lamindb/pull/495) | [falexwolf](https://github.com/falexwolf) | 2023-02-16 | 0.27.1
🚸 Add relationship between `DFolder` and `DObject` | [494](https://github.com/laminlabs/lamindb/pull/494) | [bpenteado](https://github.com/bpenteado) | 2023-02-16 |
✨ Confirm dialog for deleting data from storage | [493](https://github.com/laminlabs/lamindb/pull/493) | [sunnyosun](https://github.com/sunnyosun) | 2023-02-16 |
✨ Ingest existing data from configured local storage | [491](https://github.com/laminlabs/lamindb/pull/491) | [sunnyosun](https://github.com/sunnyosun) | 2023-02-16 |
🚸 Proper client server check for `lndb` | [492](https://github.com/laminlabs/lamindb/pull/492) | [falexwolf](https://github.com/falexwolf) | 2023-02-16 |
⬆️  Upgrade `nbproject` to 0.8.2 | [490](https://github.com/laminlabs/lamindb/pull/490) | [bpenteado](https://github.com/bpenteado) | 2023-02-16 |
⬆️  Upgrade and rename `lndb_setup` to `lndb` (v0.32.4) <span class="badge badge-warning">Breaking</span> | [487](https://github.com/laminlabs/lamindb/pull/487) | [bpenteado](https://github.com/bpenteado) | 2023-02-14 | 0.27.0
🚑 Fix tracking of added records during `ln.save()` | [489](https://github.com/laminlabs/lamindb/pull/489) | [bpenteado](https://github.com/bpenteado) | 2023-02-13 |
✨ Added `is_run_input` param to `DObject.load()` | [488](https://github.com/laminlabs/lamindb/pull/488) | [sunnyosun](https://github.com/sunnyosun) | 2023-02-13 |
🎨 Added `ln.settings.track_run_inputs_upon_load` | [486](https://github.com/laminlabs/lamindb/pull/486) | [sunnyosun](https://github.com/sunnyosun) | 2023-02-13 |
🎨 Added zarr tests back and cleaned up faq | [485](https://github.com/laminlabs/lamindb/pull/485) | [sunnyosun](https://github.com/sunnyosun) | 2023-02-13 |
🔥 Drop populating runin and tracking usage upon load | [484](https://github.com/laminlabs/lamindb/pull/484) | [sunnyosun](https://github.com/sunnyosun) | 2023-02-13 |
🚸 Better `ln.nb.header()` auto-retrieval error message | [483](https://github.com/laminlabs/lamindb/pull/483) | [falexwolf](https://github.com/falexwolf) | 2023-02-09 |
🚸  Make DObject upload ACID | [476](https://github.com/laminlabs/lamindb/pull/476) | [bpenteado](https://github.com/bpenteado) | 2023-02-08 |
📝 Add notebook on Registry lazy loading behavior to FAQ | [472](https://github.com/laminlabs/lamindb/pull/472) | [bpenteado](https://github.com/bpenteado) | 2023-02-06 |
🎨 Robust generation of `DObject._filekey` | [481](https://github.com/laminlabs/lamindb/pull/481) | [sunnyosun](https://github.com/sunnyosun) | 2023-02-06 | 0.26.1
🎨 Added erroring behavior when file doesn't exist for `ln.delete` | [480](https://github.com/laminlabs/lamindb/pull/480) | [sunnyosun](https://github.com/sunnyosun) | 2023-02-06 |
📝 Removed extra fields in dev.datasets.pbmc68k | [479](https://github.com/laminlabs/lamindb/pull/479) | [sunnyosun](https://github.com/sunnyosun) | 2023-02-06 |
➖ Removed lnbfx and fix CI | [478](https://github.com/laminlabs/lamindb/pull/478) | [sunnyosun](https://github.com/sunnyosun) | 2023-02-06 |
📝 Query book | [470](https://github.com/laminlabs/lamindb/pull/470) | [sunnyosun](https://github.com/sunnyosun) | 2023-02-06 |
🩹 Print dobject name for zarr upload | [475](https://github.com/laminlabs/lamindb/pull/475) | [Koncopd](https://github.com/Koncopd) | 2023-02-02 |
🐛 Fix load | [474](https://github.com/laminlabs/lamindb/pull/474) | [Koncopd](https://github.com/Koncopd) | 2023-02-02 |
🔥 Disable Registry relationship preview | [473](https://github.com/laminlabs/lamindb/pull/473) | [bpenteado](https://github.com/bpenteado) | 2023-02-02 |
✨ Allow ingesting existing data in the cloud | [471](https://github.com/laminlabs/lamindb/pull/471) | [sunnyosun](https://github.com/sunnyosun) | 2023-02-02 |
🐛 Correct filepath in header | [469](https://github.com/laminlabs/lamindb/pull/469) | [Koncopd](https://github.com/Koncopd) | 2023-02-01 |
🚸 Add post-setup settings manager, error on duplicate insert <span class="badge badge-warning">Breaking</span> | [466](https://github.com/laminlabs/lamindb/pull/466) | [falexwolf](https://github.com/falexwolf) | 2023-02-01 | 0.26.0
🐛 Fix fallback for notebook name | [463](https://github.com/laminlabs/lamindb/pull/463) | [falexwolf](https://github.com/falexwolf) | 2023-01-30 | 0.25.7
🔥 Remove `lns.DObject` <span class="badge badge-warning">Breaking</span> | [462](https://github.com/laminlabs/lamindb/pull/462) | [sunnyosun](https://github.com/sunnyosun) | 2023-01-30 |
⬆️ Upgrade lndb-setup 0.30.11 | [461](https://github.com/laminlabs/lamindb/pull/461) | [sunnyosun](https://github.com/sunnyosun) | 2023-01-30 | 0.25.6
⬆️ Cleaned up dependencies so that it's not as redundant to lndb-setup | [460](https://github.com/laminlabs/lamindb/pull/460) | [sunnyosun](https://github.com/sunnyosun) | 2023-01-27 |
➕ Bring back lnfbx | [459](https://github.com/laminlabs/lamindb/pull/459) | [falexwolf](https://github.com/falexwolf) | 2023-01-26 | 0.25.5
📝 Prettier section headings in the docs | [456](https://github.com/laminlabs/lamindb/pull/456) | [sunnyosun](https://github.com/sunnyosun) | 2023-01-26 |
📌 Pin `s3fs` and `gcsfs` to the latest versions | [458](https://github.com/laminlabs/lamindb/pull/458) | [falexwolf](https://github.com/falexwolf) | 2023-01-26 | 0.25.4
🚸 Auto-populate relationship-associated foreign key fields | [457](https://github.com/laminlabs/lamindb/pull/457) | [bpenteado](https://github.com/bpenteado) | 2023-01-26 | 0.25.3
🐛  Fix strict type checking for relationships | [455](https://github.com/laminlabs/lamindb/pull/455) | [bpenteado](https://github.com/bpenteado) | 2023-01-24 | 0.25.2
🩺 Increase migrations testing robustness postgres | [454](https://github.com/laminlabs/lamindb/pull/454) | [falexwolf](https://github.com/falexwolf) | 2023-01-24 |
📝 Remove linked-select notebook | [453](https://github.com/laminlabs/lamindb/pull/453) | [falexwolf](https://github.com/falexwolf) | 2023-01-24 |
⬆️ Upgrade lnschema-core to 0.25.1 | [452](https://github.com/laminlabs/lamindb/pull/452) | [bpenteado](https://github.com/bpenteado) | 2023-01-23 | 0.25.1
✨ Add explicit remote sqlite instance locking to write operations | [447](https://github.com/laminlabs/lamindb/pull/447) | [Koncopd](https://github.com/Koncopd) | 2023-01-23 |
♻️ Refactored FAQ | [448](https://github.com/laminlabs/lamindb/pull/448) | [sunnyosun](https://github.com/sunnyosun) | 2023-01-23 |
📝 Create data validation FAQ | [451](https://github.com/laminlabs/lamindb/pull/451) | [bpenteado](https://github.com/bpenteado) | 2023-01-23 |
➖ Remove s3fs dependency | [450](https://github.com/laminlabs/lamindb/pull/450) | [fredericenard](https://github.com/fredericenard) | 2023-01-23 |
⬆️ Upgrade lndb-setup to 0.30.8 | [449](https://github.com/laminlabs/lamindb/pull/449) | [fredericenard](https://github.com/fredericenard) | 2023-01-23 |
🩹 Better treat edge cases upon signup, login, failed instance loading | [446](https://github.com/laminlabs/lamindb/pull/446) | [falexwolf](https://github.com/falexwolf) | 2023-01-20 | 0.25.0
🚸 Introduce data validation on the Registry level | [445](https://github.com/laminlabs/lamindb/pull/445) | [bpenteado](https://github.com/bpenteado) | 2023-01-20 |
♻️ Reorganize quickstart and get-started | [444](https://github.com/laminlabs/lamindb/pull/444) | [sunnyosun](https://github.com/sunnyosun) | 2023-01-20 |
📝 Refactor init guide and show bionty versions in guide | [443](https://github.com/laminlabs/lamindb/pull/443) | [sunnyosun](https://github.com/sunnyosun) | 2023-01-18 | 0.24.6
⬆️ Upgrade wetlab schema to 0.13.3 | [442](https://github.com/laminlabs/lamindb/pull/442) | [sunnyosun](https://github.com/sunnyosun) | 2023-01-17 | 0.24.5
🐛 Fix taxon_id type and upgrade bionty | [441](https://github.com/laminlabs/lamindb/pull/441) | [sunnyosun](https://github.com/sunnyosun) | 2023-01-17 | 0.24.4
⬆️ Upgrade to lnschema-bionty 0.6.7 | [440](https://github.com/laminlabs/lamindb/pull/440) | [sunnyosun](https://github.com/sunnyosun) | 2023-01-17 |
⬆️ Upgrade to lndb-setup 0.30.6 | [439](https://github.com/laminlabs/lamindb/pull/439) | [sunnyosun](https://github.com/sunnyosun) | 2023-01-16 | 0.24.3
⬆️ Upgrade to lndb_setup==0.30.5 | [438](https://github.com/laminlabs/lamindb/pull/438) | [sunnyosun](https://github.com/sunnyosun) | 2023-01-16 | 0.24.2
⬆️ Upgrade lndb-setup to 0.30.4 | [437](https://github.com/laminlabs/lamindb/pull/437) | [sunnyosun](https://github.com/sunnyosun) | 2023-01-16 | 0.24.1
⬆️ Upgrade to lndb-setup 0.30.2 | [436](https://github.com/laminlabs/lamindb/pull/436) | [falexwolf](https://github.com/falexwolf) | 2023-01-16 | 0.24.0
🚸 Better hash exception | [434](https://github.com/laminlabs/lamindb/pull/434) | [falexwolf](https://github.com/falexwolf) | 2023-01-12 |
🚸 Safer session behavior 2/2 <span class="badge badge-warning">Breaking</span> | [432](https://github.com/laminlabs/lamindb/pull/432) | [falexwolf](https://github.com/falexwolf) | 2023-01-12 | 0.23.0
👷 Extend CI to py3.8-3.10 | [431](https://github.com/laminlabs/lamindb/pull/431) | [sunnyosun](https://github.com/sunnyosun) | 2023-01-12 | 0.22.5
🚸 Safer session behavior 1/2 | [430](https://github.com/laminlabs/lamindb/pull/430) | [falexwolf](https://github.com/falexwolf) | 2023-01-11 |
⬆️ Upgrade lndb-setup | [428](https://github.com/laminlabs/lamindb/pull/428) | [fredericenard](https://github.com/fredericenard) | 2023-01-10 |
📝 Re-arrange notebooks | [427](https://github.com/laminlabs/lamindb/pull/427) | [falexwolf](https://github.com/falexwolf) | 2023-01-09 |
🚸 Make `ln.nb.header()` more robust | [426](https://github.com/laminlabs/lamindb/pull/426) | [falexwolf](https://github.com/falexwolf) | 2023-01-09 | 0.22.4
📝 Improving wording of definitions | [424](https://github.com/laminlabs/lamindb/pull/424) | [falexwolf](https://github.com/falexwolf) | 2023-01-08 |
📝 Fixes for lndocs upgrade | [423](https://github.com/laminlabs/lamindb/pull/423) | [falexwolf](https://github.com/falexwolf) | 2023-01-08 |
♻️ Refactor guide | [422](https://github.com/laminlabs/lamindb/pull/422) | [sunnyosun](https://github.com/sunnyosun) | 2023-01-08 |
⬆️ Upgrade `lnschema-bionty` to 0.6.5 | [421](https://github.com/laminlabs/lamindb/pull/421) | [falexwolf](https://github.com/falexwolf) | 2023-01-05 | 0.22.3
⬆️ Upgrade to `lndb-setup` 0.28.1 | [420](https://github.com/laminlabs/lamindb/pull/420) | [falexwolf](https://github.com/falexwolf) | 2023-01-05 |
Fix typos | [418](https://github.com/laminlabs/lamindb/pull/418) | [Zethson](https://github.com/Zethson) | 2023-01-02 |
🐛 Fix bugs in `lndb set` & `lndb info` | [415](https://github.com/laminlabs/lamindb/pull/415) | [falexwolf](https://github.com/falexwolf) | 2022-12-22 | 0.22.2
✅ Error behavior for ingest | [414](https://github.com/laminlabs/lamindb/pull/414) | [falexwolf](https://github.com/falexwolf) | 2022-12-20 |
🎨 Simplify | [413](https://github.com/laminlabs/lamindb/pull/413) | [sunnyosun](https://github.com/sunnyosun) | 2022-12-16 |
⬆️ Upgrade lndb-setup | [412](https://github.com/laminlabs/lamindb/pull/412) | [falexwolf](https://github.com/falexwolf) | 2022-12-16 | 0.22.1
🚸 Better CLI & logging | [411](https://github.com/laminlabs/lamindb/pull/411) | [falexwolf](https://github.com/falexwolf) | 2022-12-15 | 0.22.0
⬆️ Updated setup | [409](https://github.com/laminlabs/lamindb/pull/409) | [sunnyosun](https://github.com/sunnyosun) | 2022-12-15 | 0.21.5
🚸 Make `nb.run` and `nb.notebook` public | [408](https://github.com/laminlabs/lamindb/pull/408) | [falexwolf](https://github.com/falexwolf) | 2022-12-14 |
⬆️ Upgrade wetlab | [407](https://github.com/laminlabs/lamindb/pull/407) | [bpenteado](https://github.com/bpenteado) | 2022-12-13 | 0.21.4
⬆️ Upgrade lndb-setup | [405](https://github.com/laminlabs/lamindb/pull/405) | [falexwolf](https://github.com/falexwolf) | 2022-12-13 | 0.21.3
⬆️ Upgrade lndb-setup | [404](https://github.com/laminlabs/lamindb/pull/404) | [sunnyosun](https://github.com/sunnyosun) | 2022-12-13 | 0.21.2
⬆️ Upgrade lndb-setup | [403](https://github.com/laminlabs/lamindb/pull/403) | [falexwolf](https://github.com/falexwolf) | 2022-12-12 | 0.21.1
📝 Fix docs | [402](https://github.com/laminlabs/lamindb/pull/402) | [falexwolf](https://github.com/falexwolf) | 2022-12-09 |
🎨 Integrate `ln.record` into `lns.DObject` → `ln.DObject` | [400](https://github.com/laminlabs/lamindb/pull/400) | [sunnyosun](https://github.com/sunnyosun) | 2022-12-09 | 0.21.0
⬆️ Upgrade setup, core & wetlab schema | [398](https://github.com/laminlabs/lamindb/pull/398) | [falexwolf](https://github.com/falexwolf) | 2022-12-08 | 0.20.0
🔥 Drop all logic related to dynamic settings | [397](https://github.com/laminlabs/lamindb/pull/397) | [fredericenard](https://github.com/fredericenard) | 2022-12-08 |
⬆️ Upgrade wetlab | [395](https://github.com/laminlabs/lamindb/pull/395) | [bpenteado](https://github.com/bpenteado) | 2022-12-06 | 0.19.4
⬆️ Upgrade wetlab | [394](https://github.com/laminlabs/lamindb/pull/394) | [bpenteado](https://github.com/bpenteado) | 2022-12-06 | 0.19.3
⬆️ Upgrade wetlab | [393](https://github.com/laminlabs/lamindb/pull/393) | [sunnyosun](https://github.com/sunnyosun) | 2022-12-06 | 0.19.2
🐛 Fix view | [392](https://github.com/laminlabs/lamindb/pull/392) | [falexwolf](https://github.com/falexwolf) | 2022-12-05 | 0.19.1
🎨 Enable inheriting wetlab schemas | [391](https://github.com/laminlabs/lamindb/pull/391) | [falexwolf](https://github.com/falexwolf) | 2022-12-05 | 0.19.0
✅ Better tests for features hashing | [390](https://github.com/laminlabs/lamindb/pull/390) | [falexwolf](https://github.com/falexwolf) | 2022-12-04 |
⬆️ Upgrade lndb-setup | [389](https://github.com/laminlabs/lamindb/pull/389) | [fredericenard](https://github.com/fredericenard) | 2022-12-04 | 0.18.9
✨ Check duplication before inserting records | [387](https://github.com/laminlabs/lamindb/pull/387) | [sunnyosun](https://github.com/sunnyosun) | 2022-11-30 | 0.18.8
🚸 Do not autoflush upon select | [386](https://github.com/laminlabs/lamindb/pull/386) | [falexwolf](https://github.com/falexwolf) | 2022-11-30 | 0.18.7
🐛 Fix bug in schema module name lookup | [384](https://github.com/laminlabs/lamindb/pull/384) | [falexwolf](https://github.com/falexwolf) | 2022-11-29 |
⬆️ Updated wetlab | [383](https://github.com/laminlabs/lamindb/pull/383) | [sunnyosun](https://github.com/sunnyosun) | 2022-11-29 | 0.18.6
⬆️ Upgrade `lnschema-core` and `lndb-setup` | [382](https://github.com/laminlabs/lamindb/pull/382) | [falexwolf](https://github.com/falexwolf) | 2022-11-28 | 0.18.5
✨ Improve lazy selectors | [375](https://github.com/laminlabs/lamindb/pull/375) | [Koncopd](https://github.com/Koncopd) | 2022-11-28 |
⬆️ Update wetlab schema | [381](https://github.com/laminlabs/lamindb/pull/381) | [sunnyosun](https://github.com/sunnyosun) | 2022-11-28 | 0.18.4
⬆️ Upgrade wetlab schema | [380](https://github.com/laminlabs/lamindb/pull/380) | [sunnyosun](https://github.com/sunnyosun) | 2022-11-28 | 0.18.3
⬆️ Upgrade `lndb-setup` & `lnschema-core` | [379](https://github.com/laminlabs/lamindb/pull/379) | [fredericenard](https://github.com/fredericenard) | 2022-11-28 | 0.18.2
⬆️ Upgrade lndb-setup to 0.18.1 | [378](https://github.com/laminlabs/lamindb/pull/378) | [falexwolf](https://github.com/falexwolf) | 2022-11-25 |
🔥 Removed bioreadout | [373](https://github.com/laminlabs/lamindb/pull/373) | [sunnyosun](https://github.com/sunnyosun) | 2022-11-24 | 0.18.1
✨ Add lazy selectors to `ln.subset` | [370](https://github.com/laminlabs/lamindb/pull/370) | [Koncopd](https://github.com/Koncopd) | 2022-11-23 |
🏗️ Persist the session | [372](https://github.com/laminlabs/lamindb/pull/372) | [falexwolf](https://github.com/falexwolf) | 2022-11-23 |
⬆️ Upgrade lndb-setup | [371](https://github.com/laminlabs/lamindb/pull/371) | [fredericenard](https://github.com/fredericenard) | 2022-11-23 |
✨ Add subset function for dobjects | [368](https://github.com/laminlabs/lamindb/pull/368) | [Koncopd](https://github.com/Koncopd) | 2022-11-23 | 0.18.0
🎨 Drop `Biometa` | [369](https://github.com/laminlabs/lamindb/pull/369) | [falexwolf](https://github.com/falexwolf) | 2022-11-22 |
🎨 Simplify schema module handling | [367](https://github.com/laminlabs/lamindb/pull/367) | [falexwolf](https://github.com/falexwolf) | 2022-11-22 |
📝 Bring back guide to linking sample-level metadata | [365](https://github.com/laminlabs/lamindb/pull/365) | [falexwolf](https://github.com/falexwolf) | 2022-11-21 | 0.17.0
📝 Add flow example back to main guide | [363](https://github.com/laminlabs/lamindb/pull/363) | [falexwolf](https://github.com/falexwolf) | 2022-11-21 |
📝 Update schema | [362](https://github.com/laminlabs/lamindb/pull/362) | [falexwolf](https://github.com/falexwolf) | 2022-11-21 |
🩹 Restore default `fsspec` for upload | [361](https://github.com/laminlabs/lamindb/pull/361) | [Koncopd](https://github.com/Koncopd) | 2022-11-20 |
🎨 Prettify API | [359](https://github.com/laminlabs/lamindb/pull/359) | [falexwolf](https://github.com/falexwolf) | 2022-11-18 |
🐛 Also return existing features | [358](https://github.com/laminlabs/lamindb/pull/358) | [falexwolf](https://github.com/falexwolf) | 2022-11-18 |
📝 Fix select gene doc | [357](https://github.com/laminlabs/lamindb/pull/357) | [falexwolf](https://github.com/falexwolf) | 2022-11-18 |
🏗️ Refactor ingest | [356](https://github.com/laminlabs/lamindb/pull/356) | [falexwolf](https://github.com/falexwolf) | 2022-11-17 |
✨ Knowledge guide | [353](https://github.com/laminlabs/lamindb/pull/353) | [sunnyosun](https://github.com/sunnyosun) | 2022-11-16 |
🎨 Separate `nb.publish` from `ingest.commit` | [355](https://github.com/laminlabs/lamindb/pull/355) | [falexwolf](https://github.com/falexwolf) | 2022-11-16 | 0.16.0
🐛 Fix data source | [354](https://github.com/laminlabs/lamindb/pull/354) | [falexwolf](https://github.com/falexwolf) | 2022-11-15 |
🚚 Move storage key to core schema | [352](https://github.com/laminlabs/lamindb/pull/352) | [falexwolf](https://github.com/falexwolf) | 2022-11-14 | 0.15.0
🐛 Fixed `species_id` in bio entity tables | [351](https://github.com/laminlabs/lamindb/pull/351) | [sunnyosun](https://github.com/sunnyosun) | 2022-11-14 |
✨ Initialize `Jupynb` `run` upon `nb.header()` | [350](https://github.com/laminlabs/lamindb/pull/350) | [falexwolf](https://github.com/falexwolf) | 2022-11-12 | 0.14.0
🏗️ Aggregate `Run` and `DTransform` | [349](https://github.com/laminlabs/lamindb/pull/349) | [falexwolf](https://github.com/falexwolf) | 2022-11-12 |
⬆️ Updated bionty | [348](https://github.com/laminlabs/lamindb/pull/348) | [sunnyosun](https://github.com/sunnyosun) | 2022-11-11 |
📝 Make notebook for `link_features` | [347](https://github.com/laminlabs/lamindb/pull/347) | [falexwolf](https://github.com/falexwolf) | 2022-11-11 | 0.13.0
✨ Join dataframes | [345](https://github.com/laminlabs/lamindb/pull/345) | [falexwolf](https://github.com/falexwolf) | 2022-11-10 |
✨ Fields in join | [344](https://github.com/laminlabs/lamindb/pull/344) | [falexwolf](https://github.com/falexwolf) | 2022-11-10 |
✅ Add tests for zarr ingest and load | [342](https://github.com/laminlabs/lamindb/pull/342) | [Koncopd](https://github.com/Koncopd) | 2022-11-10 |
🚚 Rename `PipelineRun` to `Run` | [343](https://github.com/laminlabs/lamindb/pull/343) | [falexwolf](https://github.com/falexwolf) | 2022-11-10 |
🐛 Fix `view` | [341](https://github.com/laminlabs/lamindb/pull/341) | [falexwolf](https://github.com/falexwolf) | 2022-11-09 | 0.12.1
🚚 Rename view arg | [340](https://github.com/laminlabs/lamindb/pull/340) | [falexwolf](https://github.com/falexwolf) | 2022-11-09 | 0.12.0
📝 Refactor select notebooks | [339](https://github.com/laminlabs/lamindb/pull/339) | [falexwolf](https://github.com/falexwolf) | 2022-11-08 |
📝 Simplify arg in link | [338](https://github.com/laminlabs/lamindb/pull/338) | [falexwolf](https://github.com/falexwolf) | 2022-11-08 |
♻️ Refactor select | [337](https://github.com/laminlabs/lamindb/pull/337) | [falexwolf](https://github.com/falexwolf) | 2022-11-08 |
⬆️ Update to `bionty==0.5.3` | [334](https://github.com/laminlabs/lamindb/pull/334) | [sunnyosun](https://github.com/sunnyosun) | 2022-11-06 |
🚚 Rename entries to records in link | [333](https://github.com/laminlabs/lamindb/pull/333) | [sunnyosun](https://github.com/sunnyosun) | 2022-11-04 |
📝 Improve guide | [332](https://github.com/laminlabs/lamindb/pull/332) | [falexwolf](https://github.com/falexwolf) | 2022-11-04 |
🎨 Migrate to the new schema modules | [331](https://github.com/laminlabs/lamindb/pull/331) | [falexwolf](https://github.com/falexwolf) | 2022-11-04 | 0.11.0
🔥 Removed `link_biometa` | [330](https://github.com/laminlabs/lamindb/pull/330) | [sunnyosun](https://github.com/sunnyosun) | 2022-11-03 |
🎨 Removed `get`, re-organized API docs | [329](https://github.com/laminlabs/lamindb/pull/329) | [sunnyosun](https://github.com/sunnyosun) | 2022-11-03 |
🚚 Move `.db` API to root level | [328](https://github.com/laminlabs/lamindb/pull/328) | [sunnyosun](https://github.com/sunnyosun) | 2022-11-03 |
🚚 Remove A prefix from notebooks | [327](https://github.com/laminlabs/lamindb/pull/327) | [falexwolf](https://github.com/falexwolf) | 2022-11-03 |
✨ Storage related features | [322](https://github.com/laminlabs/lamindb/pull/322) | [Koncopd](https://github.com/Koncopd) | 2022-11-03 |
🎨 Simplify generating records | [326](https://github.com/laminlabs/lamindb/pull/326) | [sunnyosun](https://github.com/sunnyosun) | 2022-11-02 |
🐛 Fixed bug in linking features | [325](https://github.com/laminlabs/lamindb/pull/325) | [sunnyosun](https://github.com/sunnyosun) | 2022-11-01 |
🍱 Added `anndata_mouse_sc_lymph_node` | [324](https://github.com/laminlabs/lamindb/pull/324) | [sunnyosun](https://github.com/sunnyosun) | 2022-11-01 |
🚚 Rename `schema._table.Table` to `table_meta` | [323](https://github.com/laminlabs/lamindb/pull/323) | [sunnyosun](https://github.com/sunnyosun) | 2022-11-01 |
🎨 Added LinkFeatureToKnowledgeTable | [320](https://github.com/laminlabs/lamindb/pull/320) | [sunnyosun](https://github.com/sunnyosun) | 2022-10-31 |
⬆️ Update lndb setup | [321](https://github.com/laminlabs/lamindb/pull/321) | [fredericenard](https://github.com/fredericenard) | 2022-10-27 | 0.10.0
⬆️ Pinned bionty version | [319](https://github.com/laminlabs/lamindb/pull/319) | [sunnyosun](https://github.com/sunnyosun) | 2022-10-26 |
⬆️ Updated to lnschema_bionty 0.4.4 | [317](https://github.com/laminlabs/lamindb/pull/317) | [sunnyosun](https://github.com/sunnyosun) | 2022-10-24 |
🔥 Remove lndb-hub import | [318](https://github.com/laminlabs/lamindb/pull/318) | [fredericenard](https://github.com/fredericenard) | 2022-10-24 |
✨ Enable to get metadata from `InstanceSettingsStore` | [310](https://github.com/laminlabs/lamindb/pull/310) | [fredericenard](https://github.com/fredericenard) | 2022-10-24 |
📝 Rename row to record when not yet added to the DB | [316](https://github.com/laminlabs/lamindb/pull/316) | [falexwolf](https://github.com/falexwolf) | 2022-10-23 |
⬆️ Update lndb hub version | [315](https://github.com/laminlabs/lamindb/pull/315) | [fredericenard](https://github.com/fredericenard) | 2022-10-23 |
🎨 Replace `insert` and `update` with `add` | [308](https://github.com/laminlabs/lamindb/pull/308) | [falexwolf](https://github.com/falexwolf) | 2022-10-22 |
⬆️ Upgrade `lnbfx` to 0.4.5 | [311](https://github.com/laminlabs/lamindb/pull/311) | [bpenteado](https://github.com/bpenteado) | 2022-10-22 |
📝 Update postgres faq notebook | [314](https://github.com/laminlabs/lamindb/pull/314) | [bpenteado](https://github.com/bpenteado) | 2022-10-22 |
🚚 Moved the rds notebooks to `rnd-demo` repo | [309](https://github.com/laminlabs/lamindb/pull/309) | [sunnyosun](https://github.com/sunnyosun) | 2022-10-21 |
✨ Added `ln.link` to populate link tables given two table entries | [307](https://github.com/laminlabs/lamindb/pull/307) | [sunnyosun](https://github.com/sunnyosun) | 2022-10-21 |
🔥 Remove hub import | [301](https://github.com/laminlabs/lamindb/pull/301) | [fredericenard](https://github.com/fredericenard) | 2022-10-21 |
🎨 Overhaul `select` and add `get` | [300](https://github.com/laminlabs/lamindb/pull/300) | [falexwolf](https://github.com/falexwolf) | 2022-10-21 |
🩹 Skip nc_evolutions table created by nocodb | [302](https://github.com/laminlabs/lamindb/pull/302) | [sunnyosun](https://github.com/sunnyosun) | 2022-10-21 |
✨ Added `knowledge` module | [299](https://github.com/laminlabs/lamindb/pull/299) | [sunnyosun](https://github.com/sunnyosun) | 2022-10-20 |
⬆️ Upgrade wetlab schema | [297](https://github.com/laminlabs/lamindb/pull/297) | [sunnyosun](https://github.com/sunnyosun) | 2022-10-20 | 0.9.6
⬆️ Added `dset` and `project` tables to core | [296](https://github.com/laminlabs/lamindb/pull/296) | [sunnyosun](https://github.com/sunnyosun) | 2022-10-19 | 0.9.5
🩹 Fixed rds nbs | [295](https://github.com/laminlabs/lamindb/pull/295) | [sunnyosun](https://github.com/sunnyosun) | 2022-10-18 |
📝 Add prisma examples and implementation ideas | [294](https://github.com/laminlabs/lamindb/pull/294) | [fredericenard](https://github.com/fredericenard) | 2022-10-18 |
♻️ Refactor linked select | [286](https://github.com/laminlabs/lamindb/pull/286) | [sunnyosun](https://github.com/sunnyosun) | 2022-10-18 | 0.9.4
✏️ Fix typo in link table entry fetching | [292](https://github.com/laminlabs/lamindb/pull/292) | [bpenteado](https://github.com/bpenteado) | 2022-10-17 |
🐛 Fixed link via link tables | [291](https://github.com/laminlabs/lamindb/pull/291) | [sunnyosun](https://github.com/sunnyosun) | 2022-10-15 | 0.9.3
✨ Add an option to use `fsspec` for upload | [288](https://github.com/laminlabs/lamindb/pull/288) | [Koncopd](https://github.com/Koncopd) | 2022-10-15 |
✨ Load returns filepath if no in-memory format is found | [287](https://github.com/laminlabs/lamindb/pull/287) | [sunnyosun](https://github.com/sunnyosun) | 2022-10-13 | 0.9.2
🎨 Clean up dtransform_in | [commit](https://github.com/laminlabs/lamindb/commit/429abf1c12b2e518e6d65329e4c4067e96d28fec) | [falexwolf](https://github.com/falexwolf) | 2022-10-13 | 0.9.1
🎨 Continue overhaul | [285](https://github.com/laminlabs/lamindb/pull/285) | [falexwolf](https://github.com/falexwolf) | 2022-10-13 | 0.9.0
🎨 Refactor linking dobjects & base64 encode checksum | [283](https://github.com/laminlabs/lamindb/pull/283) | [falexwolf](https://github.com/falexwolf) | 2022-10-13 |
🐛 View only prints existing tables | [284](https://github.com/laminlabs/lamindb/pull/284) | [sunnyosun](https://github.com/sunnyosun) | 2022-10-12 |
✨ Added `db.view`, rename `schema.draw` to `schema.view` | [282](https://github.com/laminlabs/lamindb/pull/282) | [sunnyosun](https://github.com/sunnyosun) | 2022-10-12 | 0.8.3
⬆️ Added bioreadout lookup to guide | [281](https://github.com/laminlabs/lamindb/pull/281) | [sunnyosun](https://github.com/sunnyosun) | 2022-10-12 |
🚚 Rename query to select | [280](https://github.com/laminlabs/lamindb/pull/280) | [falexwolf](https://github.com/falexwolf) | 2022-10-12 |
🩹 Fix out of sync db warning | [279](https://github.com/laminlabs/lamindb/pull/279) | [Koncopd](https://github.com/Koncopd) | 2022-10-12 |
✨ Add streaming zarr write and streaming h5ad and zarr read | [277](https://github.com/laminlabs/lamindb/pull/277) | [Koncopd](https://github.com/Koncopd) | 2022-10-12 |
💄 Cosmetics | [276](https://github.com/laminlabs/lamindb/pull/276) | [falexwolf](https://github.com/falexwolf) | 2022-10-12 |
♻️ Refactor ingest & insert | [273](https://github.com/laminlabs/lamindb/pull/273) | [falexwolf](https://github.com/falexwolf) | 2022-10-11 |
✨ Compute checksum during ingest | [274](https://github.com/laminlabs/lamindb/pull/274) | [fredericenard](https://github.com/fredericenard) | 2022-10-11 |
🚚 Moved bfx ingestion to faq | [272](https://github.com/laminlabs/lamindb/pull/272) | [sunnyosun](https://github.com/sunnyosun) | 2022-10-11 |
⬆️ Upgrade to `lndb_setup` 0.12.0 | [271](https://github.com/laminlabs/lamindb/pull/271) | [falexwolf](https://github.com/falexwolf) | 2022-10-10 | 0.8.2
📝 Improved ingest guides, fixed linked entry bugs | [270](https://github.com/laminlabs/lamindb/pull/270) | [sunnyosun](https://github.com/sunnyosun) | 2022-10-10 |
🐛 Fixed dobject_biometa entry insertion | [269](https://github.com/laminlabs/lamindb/pull/269) | [sunnyosun](https://github.com/sunnyosun) | 2022-10-10 | 0.8.1
📝 Update test users in docs | [268](https://github.com/laminlabs/lamindb/pull/268) | [falexwolf](https://github.com/falexwolf) | 2022-10-10 |
⬆️ Upgrade `lndb_setup` to 0.11.0 and `nbproject` to 0.7.0 | [267](https://github.com/laminlabs/lamindb/pull/267) | [falexwolf](https://github.com/falexwolf) | 2022-10-10 | 0.8.0
🎨 `insert.from_list` accepts entries | [266](https://github.com/laminlabs/lamindb/pull/266) | [sunnyosun](https://github.com/sunnyosun) | 2022-10-10 |
⬆️ Upgrade `nbproject` version to 0.5.5 | [265](https://github.com/laminlabs/lamindb/pull/265) | [Koncopd](https://github.com/Koncopd) | 2022-10-10 |
⬆️ Upgrade `lndb_setup` version to 0.10.1 | [264](https://github.com/laminlabs/lamindb/pull/264) | [fredericenard](https://github.com/fredericenard) | 2022-10-10 |
📝 Improved docs of `ingest`, reorganized file structure | [262](https://github.com/laminlabs/lamindb/pull/262) | [sunnyosun](https://github.com/sunnyosun) | 2022-10-10 | 0.7.2
💥 Refactor `ingest`, new pipeline ingestion logic, postgres test | [257](https://github.com/laminlabs/lamindb/pull/257) | [sunnyosun](https://github.com/sunnyosun) | 2022-10-08 | 0.7.1
📝 Overhauled get-started | [259](https://github.com/laminlabs/lamindb/pull/259) | [falexwolf](https://github.com/falexwolf) | 2022-10-07 |
🚸 Check for existence before deletion | [258](https://github.com/laminlabs/lamindb/pull/258) | [falexwolf](https://github.com/falexwolf) | 2022-10-07 |
🎨 Make ingest a static class | [256](https://github.com/laminlabs/lamindb/pull/256) | [sunnyosun](https://github.com/sunnyosun) | 2022-10-07 | 0.7.0
⬆️ Upgrade to core schema 0.10.0 | [255](https://github.com/laminlabs/lamindb/pull/255) | [falexwolf](https://github.com/falexwolf) | 2022-10-07 |
✨ New ingest API | [254](https://github.com/laminlabs/lamindb/pull/254) | [sunnyosun](https://github.com/sunnyosun) | 2022-10-05 | 0.6.0
📝 Overhaul documentation | [253](https://github.com/laminlabs/lamindb/pull/253) | [falexwolf](https://github.com/falexwolf) | 2022-10-05 |
🔧 Update lndb setup to 0.9.4 | [252](https://github.com/laminlabs/lamindb/pull/252) | [fredericenard](https://github.com/fredericenard) | 2022-10-04 |
🔧 Enable setup outside cli | [250](https://github.com/laminlabs/lamindb/pull/250) | [fredericenard](https://github.com/fredericenard) | 2022-10-04 |
⬆️ Upgrade lndb_setup and bionty | [249](https://github.com/laminlabs/lamindb/pull/249) | [falexwolf](https://github.com/falexwolf) | 2022-10-03 | 0.5.0 0.5.0
🚸 Check for migrations upon import | [247](https://github.com/laminlabs/lamindb/pull/247) | [falexwolf](https://github.com/falexwolf) | 2022-10-03 |
🔊 Raise warnings for unpopulated columns from insert.from_df | [248](https://github.com/laminlabs/lamindb/pull/248) | [sunnyosun](https://github.com/sunnyosun) | 2022-10-03 |
⬆️ Upgrade `lndb_setup` to 0.8.3 | [246](https://github.com/laminlabs/lamindb/pull/246) | [fredericenard](https://github.com/fredericenard) | 2022-10-03 |
♻️ Inherit `IngestObject` and `IngestPipelineRun` from `IngestEntity` | [243](https://github.com/laminlabs/lamindb/pull/243) | [bpenteado](https://github.com/bpenteado) | 2022-10-01 | 0.4.1
⬆️ Upgrade to lnschema_core 0.9.0 | [241](https://github.com/laminlabs/lamindb/pull/241) | [falexwolf](https://github.com/falexwolf) | 2022-09-30 |
🎨 Add type annotation to ingest | [240](https://github.com/laminlabs/lamindb/pull/240) | [sunnyosun](https://github.com/sunnyosun) | 2022-09-30 |
♻️ Generalize pipeline ingestion | [237](https://github.com/laminlabs/lamindb/pull/237) | [bpenteado](https://github.com/bpenteado) | 2022-09-30 |
✨ Added `one_or_none`, fixed api links | [238](https://github.com/laminlabs/lamindb/pull/238) | [sunnyosun](https://github.com/sunnyosun) | 2022-09-29 |
🎨 Get db metadata as a dictionary | [233](https://github.com/laminlabs/lamindb/pull/233) | [fredericenard](https://github.com/fredericenard) | 2022-09-29 |
🎨 Added `.df()` as an option to return select results | [236](https://github.com/laminlabs/lamindb/pull/236) | [sunnyosun](https://github.com/sunnyosun) | 2022-09-29 |
⬆️ Upgrade pkg versions | [234](https://github.com/laminlabs/lamindb/pull/234) | [sunnyosun](https://github.com/sunnyosun) | 2022-09-29 |
🚸 Test notebook integrity before anything else | [232](https://github.com/laminlabs/lamindb/pull/232) | [falexwolf](https://github.com/falexwolf) | 2022-09-26 |
📝 Add a quickstart & update ingest | [231](https://github.com/laminlabs/lamindb/pull/231) | [falexwolf](https://github.com/falexwolf) | 2022-09-26 |
📝 Update bfx ingestion demo | [227](https://github.com/laminlabs/lamindb/pull/227) | [bpenteado](https://github.com/bpenteado) | 2022-09-26 |
🚚 Rename insert.features to insert.featureset_from_features | [230](https://github.com/laminlabs/lamindb/pull/230) | [sunnyosun](https://github.com/sunnyosun) | 2022-09-26 |
📝 Polish documentation | [229](https://github.com/laminlabs/lamindb/pull/229) | [falexwolf](https://github.com/falexwolf) | 2022-09-26 |
✨ Populate `dobject` size | [228](https://github.com/laminlabs/lamindb/pull/228) | [falexwolf](https://github.com/falexwolf) | 2022-09-26 |
🚸 Sort select results in DataFrame by source | [226](https://github.com/laminlabs/lamindb/pull/226) | [falexwolf](https://github.com/falexwolf) | 2022-09-25 |
⬆️ Upgrade to core schema 0.7.3 | [223](https://github.com/laminlabs/lamindb/pull/223) | [falexwolf](https://github.com/falexwolf) | 2022-09-25 |
🩹 Removed logging for inserting link tables | [221](https://github.com/laminlabs/lamindb/pull/221) | [sunnyosun](https://github.com/sunnyosun) | 2022-09-24 |
✨ Ingest bfx outputs with sample metadata | [215](https://github.com/laminlabs/lamindb/pull/215) | [bpenteado](https://github.com/bpenteado) | 2022-09-24 |
👷 Improve notebook test function | [218](https://github.com/laminlabs/lamindb/pull/218) | [falexwolf](https://github.com/falexwolf) | 2022-09-23 |
🧱 Improved code infra | [216](https://github.com/laminlabs/lamindb/pull/216) | [sunnyosun](https://github.com/sunnyosun) | 2022-09-23 | 0.4.0
✨ Account for multiple `filepath` suffixes | [214](https://github.com/laminlabs/lamindb/pull/214) | [falexwolf](https://github.com/falexwolf) | 2022-09-23 |
♻️ Refactored feature model ingestion | [213](https://github.com/laminlabs/lamindb/pull/213) | [sunnyosun](https://github.com/sunnyosun) | 2022-09-23 |
♻️ Refactored insert | [211](https://github.com/laminlabs/lamindb/pull/211) | [sunnyosun](https://github.com/sunnyosun) | 2022-09-22 |
⬆️ Upgrade `lndb_setup` and `lnschema_bionty` | [212](https://github.com/laminlabs/lamindb/pull/212) | [falexwolf](https://github.com/falexwolf) | 2022-09-22 |
⬆️ Upgrade to core schema 0.7.2 | [208](https://github.com/laminlabs/lamindb/pull/208) | [falexwolf](https://github.com/falexwolf) | 2022-09-21 |
🚑 `dtransform` is either `pipeline_run` or `jupynb` | [207](https://github.com/laminlabs/lamindb/pull/207) | [sunnyosun](https://github.com/sunnyosun) | 2022-09-21 | 0.3.11
🐛 Fixed column mapping | [206](https://github.com/laminlabs/lamindb/pull/206) | [sunnyosun](https://github.com/sunnyosun) | 2022-09-19 | 0.3.10
✨ Allow batch insertion | [205](https://github.com/laminlabs/lamindb/pull/205) | [sunnyosun](https://github.com/sunnyosun) | 2022-09-19 |
♻️ Refactored linked queries | [204](https://github.com/laminlabs/lamindb/pull/204) | [sunnyosun](https://github.com/sunnyosun) | 2022-09-15 |
⬆️ Update lnbfx to 0.3.5 | [199](https://github.com/laminlabs/lamindb/pull/199) | [bpenteado](https://github.com/bpenteado) | 2022-09-14 | 0.3.9
🔊 Refactor pipeline logging | [197](https://github.com/laminlabs/lamindb/pull/197) | [bpenteado](https://github.com/bpenteado) | 2022-09-14 | 0.3.8
✨ Insert unmapped features | [198](https://github.com/laminlabs/lamindb/pull/198) | [sunnyosun](https://github.com/sunnyosun) | 2022-09-14 |
♻️ Refactored linked queries | [196](https://github.com/laminlabs/lamindb/pull/196) | [sunnyosun](https://github.com/sunnyosun) | 2022-09-13 |
⬆️ Updated wetlab schema | [195](https://github.com/laminlabs/lamindb/pull/195) | [sunnyosun](https://github.com/sunnyosun) | 2022-09-13 |
♻️ Refactor pipeline logging | [193](https://github.com/laminlabs/lamindb/pull/193) | [bpenteado](https://github.com/bpenteado) | 2022-09-13 | 0.3.7
⬆️ Updated lnbfx version | [194](https://github.com/laminlabs/lamindb/pull/194) | [sunnyosun](https://github.com/sunnyosun) | 2022-09-12 | 0.3.6
🍱 Updated cell ranger test dir | [192](https://github.com/laminlabs/lamindb/pull/192) | [sunnyosun](https://github.com/sunnyosun) | 2022-09-12 |
🐛 Fix synchronization error in load | [191](https://github.com/laminlabs/lamindb/pull/191) | [Koncopd](https://github.com/Koncopd) | 2022-09-12 |
⬆️ Updated setup version | [189](https://github.com/laminlabs/lamindb/pull/189) | [sunnyosun](https://github.com/sunnyosun) | 2022-09-08 | 0.3.5
♻️ Split ingestion logic into IngestObject and IngestPipeline | [188](https://github.com/laminlabs/lamindb/pull/188) | [bpenteado](https://github.com/bpenteado) | 2022-09-06 |
📝 Add flow data ingestion example | [186](https://github.com/laminlabs/lamindb/pull/186) | [sunnyosun](https://github.com/sunnyosun) | 2022-09-06 |
💄 Added sidebar to guide | [185](https://github.com/laminlabs/lamindb/pull/185) | [sunnyosun](https://github.com/sunnyosun) | 2022-09-06 |
♻️ Split pipeline ingestion from non-pipeline ingestion | [183](https://github.com/laminlabs/lamindb/pull/183) | [bpenteado](https://github.com/bpenteado) | 2022-09-06 |
⬆️ Updated schema module versions | [184](https://github.com/laminlabs/lamindb/pull/184) | [sunnyosun](https://github.com/sunnyosun) | 2022-09-06 |
♻️ Refactored `select` | [182](https://github.com/laminlabs/lamindb/pull/182) | [sunnyosun](https://github.com/sunnyosun) | 2022-09-05 |
🚚 Rename `return_df` to `as_df` in `select` | [181](https://github.com/laminlabs/lamindb/pull/181) | [sunnyosun](https://github.com/sunnyosun) | 2022-09-05 |
🚚 Rename `guides` to `faq` & `tutorials` to `guide` | [180](https://github.com/laminlabs/lamindb/pull/180) | [falexwolf](https://github.com/falexwolf) | 2022-09-05 |
♻️ Refactored guide | [179](https://github.com/laminlabs/lamindb/pull/179) | [sunnyosun](https://github.com/sunnyosun) | 2022-09-03 |
📝 Prettify documentation | [178](https://github.com/laminlabs/lamindb/pull/178) | [falexwolf](https://github.com/falexwolf) | 2022-09-02 |
Update lndb_hub version | [176](https://github.com/laminlabs/lamindb/pull/176) | [fredericenard](https://github.com/fredericenard) | 2022-09-01 |
✨ Allow selecting dobjects by biological entities | [175](https://github.com/laminlabs/lamindb/pull/175) | [sunnyosun](https://github.com/sunnyosun) | 2022-09-01 | 0.3.4
⬆️ Update to sqm 0.0.8 to silence the warnings | [174](https://github.com/laminlabs/lamindb/pull/174) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-30 |
🤡 R&D team simulation | [172](https://github.com/laminlabs/lamindb/pull/172) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-30 |
🚸 Fix pipeline ingestion logging | [158](https://github.com/laminlabs/lamindb/pull/158) | [bpenteado](https://github.com/bpenteado) | 2022-08-30 |
🚚 Move problems page to lamin-profile | [171](https://github.com/laminlabs/lamindb/pull/171) | [falexwolf](https://github.com/falexwolf) | 2022-08-30 |
🩹 Uncomment out sharing tests | [168](https://github.com/laminlabs/lamindb/pull/168) | [fredericenard](https://github.com/fredericenard) | 2022-08-30 |
🔥 Removed examples dir | [170](https://github.com/laminlabs/lamindb/pull/170) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-30 |
🚸 Improve delete function | [169](https://github.com/laminlabs/lamindb/pull/169) | [falexwolf](https://github.com/falexwolf) | 2022-08-30 |
🎨 Remove `track` submodule, move session, rename biogram to erdiagram | [167](https://github.com/laminlabs/lamindb/pull/167) | [falexwolf](https://github.com/falexwolf) | 2022-08-29 |
🚚 Rename `.do` to `.db` | [166](https://github.com/laminlabs/lamindb/pull/166) | [falexwolf](https://github.com/falexwolf) | 2022-08-29 |
🎨 Move header call to `nb`, re-export all of nbproject | [165](https://github.com/laminlabs/lamindb/pull/165) | [falexwolf](https://github.com/falexwolf) | 2022-08-29 |
👷 Get rid of sqm warnings | [164](https://github.com/laminlabs/lamindb/pull/164) | [falexwolf](https://github.com/falexwolf) | 2022-08-29 |
🎨 Simplify loading data | [163](https://github.com/laminlabs/lamindb/pull/163) | [falexwolf](https://github.com/falexwolf) | 2022-08-29 |
👷 Allow stripping notebooks, upgrade `nbproject_test`  | [162](https://github.com/laminlabs/lamindb/pull/162) | [falexwolf](https://github.com/falexwolf) | 2022-08-29 | 0.3.3
✨ Use `cell_marker` feature model for flow data | [161](https://github.com/laminlabs/lamindb/pull/161) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-28 |
⬆️ Updated to bionty 0.2.2 | [160](https://github.com/laminlabs/lamindb/pull/160) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-28 |
⬆️ Updated to sqm 0.0.7 | [159](https://github.com/laminlabs/lamindb/pull/159) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-28 |
⬆️ Updated dependencies | [157](https://github.com/laminlabs/lamindb/pull/157) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-28 |
⬆️ Upgrade to `lnschema_core` 0.5.1 | [156](https://github.com/laminlabs/lamindb/pull/156) | [falexwolf](https://github.com/falexwolf) | 2022-08-26 |
♻️ Update pipeline ingestion to `lnschema_core` 0.5.0 | [154](https://github.com/laminlabs/lamindb/pull/154) | [bpenteado](https://github.com/bpenteado) | 2022-08-26 |
🍱 Add scrnaseq cellranger dataset | [151](https://github.com/laminlabs/lamindb/pull/151) | [bpenteado](https://github.com/bpenteado) | 2022-08-26 |
🐛 Fix population of `dtransform` | [155](https://github.com/laminlabs/lamindb/pull/155) | [falexwolf](https://github.com/falexwolf) | 2022-08-26 |
✨ Populate `dtransform_in` | [153](https://github.com/laminlabs/lamindb/pull/153) | [falexwolf](https://github.com/falexwolf) | 2022-08-26 |
⬆️ Updated dependencies | [152](https://github.com/laminlabs/lamindb/pull/152) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-26 | 0.3.2
♻️ Refactor `lnbfx` integration | [149](https://github.com/laminlabs/lamindb/pull/149) | [bpenteado](https://github.com/bpenteado) | 2022-08-26 |
🚧 Temporary solution to extend modules | [150](https://github.com/laminlabs/lamindb/pull/150) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-25 | 0.3.1
🏗️ Make tables within schema modules configurable | [148](https://github.com/laminlabs/lamindb/pull/148) | [falexwolf](https://github.com/falexwolf) | 2022-08-25 | 0.3.0
✏️ Fixed typo in check versions | [147](https://github.com/laminlabs/lamindb/pull/147) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-24 |
🏗️ Use id to reference storage | [146](https://github.com/laminlabs/lamindb/pull/146) | [fredericenard](https://github.com/fredericenard) | 2022-08-23 |
🚚 Renamed `FeatureModel` to `LinkFeatureModel` | [145](https://github.com/laminlabs/lamindb/pull/145) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-23 |
⬆️ Upgrade to lnbfx 0.2.0 | [144](https://github.com/laminlabs/lamindb/pull/144) | [falexwolf](https://github.com/falexwolf) | 2022-08-23 |
🚚 Renamed `id` to `key` in `update` and `delete` APIs | [143](https://github.com/laminlabs/lamindb/pull/143) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-23 |
♻️ Cleaned up and added a mouse dataset | [142](https://github.com/laminlabs/lamindb/pull/142) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-23 |
⬆️ Upgrade lndb-hub to v0.5.0 | [140](https://github.com/laminlabs/lamindb/pull/140) | [fredericenard](https://github.com/fredericenard) | 2022-08-23 |
💄 Added logging to update and delete | [141](https://github.com/laminlabs/lamindb/pull/141) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-22 |
✨ Autogenerate `select`, `update` and `delete` | [139](https://github.com/laminlabs/lamindb/pull/139) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-22 |
⬆️ Upgrade core schema to v0.4.0 | [138](https://github.com/laminlabs/lamindb/pull/138) | [falexwolf](https://github.com/falexwolf) | 2022-08-22 |
✨ Track storage root location | [137](https://github.com/laminlabs/lamindb/pull/137) | [fredericenard](https://github.com/fredericenard) | 2022-08-22 |
🚚 Rename schema modules and bioinformatics module | [136](https://github.com/laminlabs/lamindb/pull/136) | [falexwolf](https://github.com/falexwolf) | 2022-08-19 |
✨ Integrate bioinformatics pipline runs | [133](https://github.com/laminlabs/lamindb/pull/133) | [bpenteado](https://github.com/bpenteado) | 2022-08-18 |
📝 Updated docs | [135](https://github.com/laminlabs/lamindb/pull/135) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-18 |
👽 Adapt to `lnschema-bionty 0.1.4` | [134](https://github.com/laminlabs/lamindb/pull/134) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-18 |
⚡ Ingest genes properly | [131](https://github.com/laminlabs/lamindb/pull/131) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-16 |
✨ New readout table, registered species | [130](https://github.com/laminlabs/lamindb/pull/130) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-12 |
🚚 Migrate `bioreader` to `bioreadout` | [129](https://github.com/laminlabs/lamindb/pull/129) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-04 |
⬆️ Upgrade to lndb_setup 0.5.0 | [128](https://github.com/laminlabs/lamindb/pull/128) | [falexwolf](https://github.com/falexwolf) | 2022-08-03 |
⬆️ Upgrade to `lnschema_core` 0.3.0 | [127](https://github.com/laminlabs/lamindb/pull/127) | [falexwolf](https://github.com/falexwolf) | 2022-08-03 |
🚚 Renamed `meta.annotate` to `db.link` | [126](https://github.com/laminlabs/lamindb/pull/126) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-03 |
👷 Cleaner CI environment switching | [125](https://github.com/laminlabs/lamindb/pull/125) | [falexwolf](https://github.com/falexwolf) | 2022-08-02 | 0.2.1
⬆️ Upgrade to `lndb_setup` 0.4.2 | [123](https://github.com/laminlabs/lamindb/pull/123) | [falexwolf](https://github.com/falexwolf) | 2022-08-02 | 0.2.0
📝 Polish the `introspect` tutorial | [122](https://github.com/laminlabs/lamindb/pull/122) | [falexwolf](https://github.com/falexwolf) | 2022-08-02 |
📝 Polish the 4 key tutorial pages | [121](https://github.com/laminlabs/lamindb/pull/121) | [falexwolf](https://github.com/falexwolf) | 2022-08-02 |
📝 Polish `get-started` & `select-data` | [119](https://github.com/laminlabs/lamindb/pull/119) | [falexwolf](https://github.com/falexwolf) | 2022-08-01 |
✨ Allow ingesting fcs files and selecting genes | [118](https://github.com/laminlabs/lamindb/pull/118) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-01 |
📝 Prettify user mentions in guide & faq | [117](https://github.com/laminlabs/lamindb/pull/117) | [falexwolf](https://github.com/falexwolf) | 2022-08-01 |
⬆️ Upgrade to `lndb_setup` 0.4.0 | [116](https://github.com/laminlabs/lamindb/pull/116) | [falexwolf](https://github.com/falexwolf) | 2022-08-01 |
🎨 Overhauled guide, renamed `load` to `select.table_as_df` | [115](https://github.com/laminlabs/lamindb/pull/115) | [sunnyosun](https://github.com/sunnyosun) | 2022-08-01 |
🩹 Some fixes of ingesting in-memory dobjects and an analysis draft | [114](https://github.com/laminlabs/lamindb/pull/114) | [falexwolf](https://github.com/falexwolf) | 2022-07-31 |
✨ Select and update metadata | [113](https://github.com/laminlabs/lamindb/pull/113) | [sunnyosun](https://github.com/sunnyosun) | 2022-07-31 |
📝 Add a problem statement | [112](https://github.com/laminlabs/lamindb/pull/112) | [falexwolf](https://github.com/falexwolf) | 2022-07-31 |
⬆️ Upgrade to `lamindb-schema` 0.3.1 | [111](https://github.com/laminlabs/lamindb/pull/111) | [falexwolf](https://github.com/falexwolf) | 2022-07-31 |
🚧 Annotate features during ingestion | [110](https://github.com/laminlabs/lamindb/pull/110) | [sunnyosun](https://github.com/sunnyosun) | 2022-07-31 |
📝 Re-write landing page | [109](https://github.com/laminlabs/lamindb/pull/109) | [falexwolf](https://github.com/falexwolf) | 2022-07-30 |
📝 Improve `get-started` and `collaborate` guide | [108](https://github.com/laminlabs/lamindb/pull/108) | [falexwolf](https://github.com/falexwolf) | 2022-07-30 |
⬆️ Upgrade to `lnschema_core` 0.2.0 | [106](https://github.com/laminlabs/lamindb/pull/106) | [falexwolf](https://github.com/falexwolf) | 2022-07-29 |
🚚 Move out `lndb_hub` | [105](https://github.com/laminlabs/lamindb/pull/105) | [falexwolf](https://github.com/falexwolf) | 2022-07-29 |
♻️ Refactor sharing on the hub | [104](https://github.com/laminlabs/lamindb/pull/104) | [falexwolf](https://github.com/falexwolf) | 2022-07-29 |
✨ Enable sharing dobjects & instances in the hub | [69](https://github.com/laminlabs/lamindb/pull/69) | [fredericenard](https://github.com/fredericenard) | 2022-07-29 |
🔥 Moved readout vocab to bioreader | [103](https://github.com/laminlabs/lamindb/pull/103) | [sunnyosun](https://github.com/sunnyosun) | 2022-07-28 |
⬆️ Upgrade `lndb_setup` to 0.3.0 | [102](https://github.com/laminlabs/lamindb/pull/102) | [falexwolf](https://github.com/falexwolf) | 2022-07-26 |
⬆️ Upgrade to lndb-cli 0.2.0 | [101](https://github.com/laminlabs/lamindb/pull/101) | [falexwolf](https://github.com/falexwolf) | 2022-07-25 |
🚚 Account for table construction in lndb-cli | [99](https://github.com/laminlabs/lamindb/pull/99) | [falexwolf](https://github.com/falexwolf) | 2022-07-25 |
🚚 Move CLI code to `lndb-cli` | [98](https://github.com/laminlabs/lamindb/pull/98) | [falexwolf](https://github.com/falexwolf) | 2022-07-24 |
🍱 Adapted to `lnschema-biology` `0.1.1` | [97](https://github.com/laminlabs/lamindb/pull/97) | [sunnyosun](https://github.com/sunnyosun) | 2022-07-24 |
🚸 Let the CLI have actual subcommands | [96](https://github.com/laminlabs/lamindb/pull/96) | [falexwolf](https://github.com/falexwolf) | 2022-07-24 |
👷 Add time out to GitHub Actions | [93](https://github.com/laminlabs/lamindb/pull/93) | [Koncopd](https://github.com/Koncopd) | 2022-07-23 |
🔥 Switched logger to use lamin-utils | [92](https://github.com/laminlabs/lamindb/pull/92) | [sunnyosun](https://github.com/sunnyosun) | 2022-07-23 |
✨  Annotate biometa when annotating genes | [90](https://github.com/laminlabs/lamindb/pull/90) | [sunnyosun](https://github.com/sunnyosun) | 2022-07-23 |
🔼 Upgrade to schema 0.2.1 | [89](https://github.com/laminlabs/lamindb/pull/89) | [falexwolf](https://github.com/falexwolf) | 2022-07-23 |
🎨 Use modular schema structure | [88](https://github.com/laminlabs/lamindb/pull/88) | [falexwolf](https://github.com/falexwolf) | 2022-07-22 |
🚚 Move db API from admin to dev | [87](https://github.com/laminlabs/lamindb/pull/87) | [falexwolf](https://github.com/falexwolf) | 2022-07-22 |
✨ Enable annotating features of dobjects | [81](https://github.com/laminlabs/lamindb/pull/81) | [sunnyosun](https://github.com/sunnyosun) | 2022-07-22 |
🚸 Offer manual way of completing a migration | [86](https://github.com/laminlabs/lamindb/pull/86) | [falexwolf](https://github.com/falexwolf) | 2022-07-22 |
🏗️ Separate settings into user vs. instance and one file per instance | [83](https://github.com/laminlabs/lamindb/pull/83) | [fredericenard](https://github.com/fredericenard) | 2022-07-22 |
🐛 Fix loading of multi-index and non-id tables | [85](https://github.com/laminlabs/lamindb/pull/85) | [falexwolf](https://github.com/falexwolf) | 2022-07-21 | 0.1.2
🚚 Rename table `interface` to `jupynb`: migrate to schema v0.1.1 | [84](https://github.com/laminlabs/lamindb/pull/84) | [falexwolf](https://github.com/falexwolf) | 2022-07-21 |
🚚 Move version check into correct __init__ | [82](https://github.com/laminlabs/lamindb/pull/82) | [falexwolf](https://github.com/falexwolf) | 2022-07-20 |
⬆️ Fix the publishing call by upgrading to nbproject 0.4.3 | [80](https://github.com/laminlabs/lamindb/pull/80) | [falexwolf](https://github.com/falexwolf) | 2022-07-19 | 0.1.1
✨ Add `schema_version` check | [79](https://github.com/laminlabs/lamindb/pull/79) | [falexwolf](https://github.com/falexwolf) | 2022-07-19 | 0.1.0
🚚 Rename `lamindb.model` to `laminln.schema` | [77](https://github.com/laminlabs/lamindb/pull/77) | [falexwolf](https://github.com/falexwolf) | 2022-07-17 |
🚚 Migrate schema out to `lamindb-schema` | [76](https://github.com/laminlabs/lamindb/pull/76) | [falexwolf](https://github.com/falexwolf) | 2022-07-17 |
✨ Version dobjects and interfaces | [75](https://github.com/laminlabs/lamindb/pull/75) | [falexwolf](https://github.com/falexwolf) | 2022-07-16 |
💄 Pretty logging | [74](https://github.com/laminlabs/lamindb/pull/74) | [sunnyosun](https://github.com/sunnyosun) | 2022-07-15 |
✨ Add `ingest.add` and `ingest.commit` | [73](https://github.com/laminlabs/lamindb/pull/73) | [sunnyosun](https://github.com/sunnyosun) | 2022-07-14 |
📝 Add an example of ingesting images | [72](https://github.com/laminlabs/lamindb/pull/72) | [sunnyosun](https://github.com/sunnyosun) | 2022-07-12 |
⬆️ Upgrade to nbproject 0.2.1 | [71](https://github.com/laminlabs/lamindb/pull/71) | [falexwolf](https://github.com/falexwolf) | 2022-07-12 |
🐛 Fix bug | [70](https://github.com/laminlabs/lamindb/pull/70) | [falexwolf](https://github.com/falexwolf) | 2022-07-12 |
📝 Added example notebook for ingesting fcs files | [66](https://github.com/laminlabs/lamindb/pull/66) | [sunnyosun](https://github.com/sunnyosun) | 2022-07-11 |
⬆️ Migrate to nbproject 0.2.0 | [68](https://github.com/laminlabs/lamindb/pull/68) | [falexwolf](https://github.com/falexwolf) | 2022-07-11 |
🚸 Auto-check integrity upon data ingestion only on Jupyter Lab | [65](https://github.com/laminlabs/lamindb/pull/65) | [falexwolf](https://github.com/falexwolf) | 2022-07-09 | 0.0.9
🚚 Rename`lndb config` to `lndb init` and rewrite get-started | [64](https://github.com/laminlabs/lamindb/pull/64) | [falexwolf](https://github.com/falexwolf) | 2022-07-09 |
🚸 Raise error upon multiple sign ups with same unconfirmed email | [63](https://github.com/laminlabs/lamindb/pull/63) | [falexwolf](https://github.com/falexwolf) | 2022-07-09 |
✨ Allow sharing instances with other users | [62](https://github.com/laminlabs/lamindb/pull/62) | [falexwolf](https://github.com/falexwolf) | 2022-07-08 |
🏗️ Improve configuration logic, flow, logging & testing | [61](https://github.com/laminlabs/lamindb/pull/61) | [falexwolf](https://github.com/falexwolf) | 2022-07-08 |
✨ Populate user metadata at sign up | [60](https://github.com/laminlabs/lamindb/pull/60) | [fredericenard](https://github.com/fredericenard) | 2022-07-07 |
✨ Unique user identity across instances | [59](https://github.com/laminlabs/lamindb/pull/59) | [falexwolf](https://github.com/falexwolf) | 2022-07-06 |
🚧 Enable storing sqlite on S3 | [58](https://github.com/laminlabs/lamindb/pull/58) | [falexwolf](https://github.com/falexwolf) | 2022-07-04 | 0.0.8
🚚 Rename table `file` to `dobject` | [57](https://github.com/laminlabs/lamindb/pull/57) | [falexwolf](https://github.com/falexwolf) | 2022-07-04 |
🚚 Rename CLI from `lamindb` to `lndb` | [56](https://github.com/laminlabs/lamindb/pull/56) | [falexwolf](https://github.com/falexwolf) | 2022-07-04 |
✨ Add `instance_name` to settings and `instance` to setup | [55](https://github.com/laminlabs/lamindb/pull/55) | [falexwolf](https://github.com/falexwolf) | 2022-07-03 | 0.0.7
♻️ Refactor file storage | [54](https://github.com/laminlabs/lamindb/pull/54) | [falexwolf](https://github.com/falexwolf) | 2022-07-03 | 0.0.6
♻️ Refactor settings | [52](https://github.com/laminlabs/lamindb/pull/52) | [falexwolf](https://github.com/falexwolf) | 2022-07-02 | 0.0.5
👷 Switch to nbproject tests & clean up logging | [51](https://github.com/laminlabs/lamindb/pull/51) | [falexwolf](https://github.com/falexwolf) | 2022-07-01 |
✨ Use nbproject publish functionality | [50](https://github.com/laminlabs/lamindb/pull/50) | [falexwolf](https://github.com/falexwolf) | 2022-07-01 |
🚚 Rename tutorial to guide | [49](https://github.com/laminlabs/lamindb/pull/49) | [falexwolf](https://github.com/falexwolf) | 2022-06-29 |
👷 Measure coverage | [48](https://github.com/laminlabs/lamindb/pull/48) | [sunnyosun](https://github.com/sunnyosun) | 2022-06-26 |
🏗️ Settings: From pydantic BaseModel to custom class | [46](https://github.com/laminlabs/lamindb/pull/46) | [falexwolf](https://github.com/falexwolf) | 2022-06-25 |
🏗️ Complete refactor of setup, settings & storage management  | [45](https://github.com/laminlabs/lamindb/pull/45) | [falexwolf](https://github.com/falexwolf) | 2022-06-25 |
⬆️ Upgrade to nbproject 0.1a3 | [43](https://github.com/laminlabs/lamindb/pull/43) | [falexwolf](https://github.com/falexwolf) | 2022-06-23 |
📝 One-page layout for data models | [42](https://github.com/laminlabs/lamindb/pull/42) | [falexwolf](https://github.com/falexwolf) | 2022-06-14 |
♻️ More explicit `_id` name for such fields | [41](https://github.com/laminlabs/lamindb/pull/41) | [falexwolf](https://github.com/falexwolf) | 2022-06-14 |
♻️ Migrate from sql to sqm everywhere | [40](https://github.com/laminlabs/lamindb/pull/40) | [falexwolf](https://github.com/falexwolf) | 2022-06-14 |
✨ Implement data access log `track.do` | [39](https://github.com/laminlabs/lamindb/pull/39) | [falexwolf](https://github.com/falexwolf) | 2022-06-12 | 0.0.4
🏗️ Re-organize API | [38](https://github.com/laminlabs/lamindb/pull/38) | [falexwolf](https://github.com/falexwolf) | 2022-06-12 |
🏗️ Set up db with sqlmodel, test int ids | [37](https://github.com/laminlabs/lamindb/pull/37) | [falexwolf](https://github.com/falexwolf) | 2022-06-11 |
🏗️ Re-designed entire API | [36](https://github.com/laminlabs/lamindb/pull/36) | [falexwolf](https://github.com/falexwolf) | 2022-06-10 |
🔥 Remove notion integration up to CLI | [35](https://github.com/laminlabs/lamindb/pull/35) | [falexwolf](https://github.com/falexwolf) | 2022-06-10 |
📝 Rename from lamindb to LaminDB & rewrite the landing page | [34](https://github.com/laminlabs/lamindb/pull/34) | [falexwolf](https://github.com/falexwolf) | 2022-06-10 | 0.0.3
🏗️ Introduce field `interface.type` | [31](https://github.com/laminlabs/lamindb/pull/31) | [falexwolf](https://github.com/falexwolf) | 2022-06-09 |
🚚  Rename global `source` table to `interface` | [30](https://github.com/laminlabs/lamindb/pull/30) | [falexwolf](https://github.com/falexwolf) | 2022-06-09 |
✨ Track python package dependencies in `source.dependency` | [27](https://github.com/laminlabs/lamindb/pull/27) | [falexwolf](https://github.com/falexwolf) | 2022-06-09 |
✨ Track title of ingesting notebook in `source.name` | [26](https://github.com/laminlabs/lamindb/pull/26) | [falexwolf](https://github.com/falexwolf) | 2022-06-09 |
🐛 Fix user & notebook ingestion, add another test dataset | [24](https://github.com/laminlabs/lamindb/pull/24) | [falexwolf](https://github.com/falexwolf) | 2022-06-08 |
✅ Fix pandas load, add tests | [23](https://github.com/laminlabs/lamindb/pull/23) | [falexwolf](https://github.com/falexwolf)
✨ Add introspection: `db.diagram()`, `db.entities()`, `db.load()` | [22](https://github.com/laminlabs/lamindb/pull/22) | [falexwolf](https://github.com/falexwolf) | |
✨ Add entity `user` | [21](https://github.com/laminlabs/lamindb/pull/21) | [falexwolf](https://github.com/falexwolf) | |
🚸 Auto-create local storage dir & cache dir | [20](https://github.com/laminlabs/lamindb/pull/20) | [falexwolf](https://github.com/falexwolf) | |
✅ Add a test for db creation & file ingestion | [18](https://github.com/laminlabs/lamindb/pull/18) | [falexwolf](https://github.com/falexwolf) | |
🏗️ Name the database file like the storage root directory | [17](https://github.com/laminlabs/lamindb/pull/17) | [falexwolf](https://github.com/falexwolf) | |
♻️ Refactor `lndb.db` | [16](https://github.com/laminlabs/lamindb/pull/16) | [falexwolf](https://github.com/falexwolf) | |
🔧 Refactor `lndb.settings` | [15](https://github.com/laminlabs/lamindb/pull/15) | [falexwolf](https://github.com/falexwolf) | |
🔥 Remove versioneer | [14](https://github.com/laminlabs/lamindb/pull/14) | [falexwolf](https://github.com/falexwolf) | |
👷 Track changes as in `cookiecutter-py` 0.3.0 | [13](https://github.com/laminlabs/lamindb/pull/13) | [falexwolf](https://github.com/falexwolf) | 2022-06-07 |
🏗️ Basic file ingestion from a notebook | [9](https://github.com/laminlabs/lamindb/pull/9) | [falexwolf](https://github.com/falexwolf) | 2022-06-06 |
🔥 Remove `sqlmodel` dependency | [8](https://github.com/laminlabs/lamindb/pull/8) | [falexwolf](https://github.com/falexwolf) | 2022-06-04 |
📄 Change format of license file | [7](https://github.com/laminlabs/lamindb/pull/7) | [falexwolf](https://github.com/falexwolf) | 2022-05-23 |
📝 Fix faq link | [6](https://github.com/laminlabs/lamindb/pull/6) | [falexwolf](https://github.com/falexwolf) | 2022-05-26 |
♻️ Update to cookiecutter 0.2.0 | [5](https://github.com/laminlabs/lamindb/pull/5) | [falexwolf](https://github.com/falexwolf) | 2022-05-23 |
💄 More narrow sidebar | [4](https://github.com/laminlabs/lamindb/pull/4) | [falexwolf](https://github.com/falexwolf) | 2022-05-11 |
💄 Switch `faq` and `api` | [3](https://github.com/laminlabs/lamindb/pull/3) | [falexwolf](https://github.com/falexwolf) | 2022-05-11 |
📝 Polish the documentation | [2](https://github.com/laminlabs/lamindb/pull/2) | [falexwolf](https://github.com/falexwolf) | 2022-05-11 |
🚚 Migrated from `lamin` repository | [1](https://github.com/laminlabs/lamindb/pull/1) | [falexwolf](https://github.com/falexwolf) | |
