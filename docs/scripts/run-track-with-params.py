import argparse
import lamindb as ln

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--input-dir", type=str)
    p.add_argument("--downsample", action="store_true")
    p.add_argument("--learning-rate", type=float)
    args = p.parse_args()
    params = {
        "input_dir": args.input_dir,
        "learning_rate": args.learning_rate,
        "downsample": args.downsample,
    }
    ln.track("JjRF4mACd9m00000", params=params)
    # your code
    ln.finish()
