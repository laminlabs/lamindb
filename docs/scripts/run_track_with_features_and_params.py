import argparse
import lamindb as ln


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--s3-folder", type=str)
    p.add_argument("--experiment", type=str)
    args = p.parse_args()
    features = {
        "s3_folder": args.s3_folder,
        "experiment": args.experiment,
    }
    ln.track(features=features, params={"example_param": 42})

    # your code

    ln.finish()
