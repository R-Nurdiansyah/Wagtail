import sys
import os

def main():
    args = sys.argv[1:]
    
    # Check if any GreenGenes-related arguments are present
    gg_keywords = ["gg2", "greengenes", "greengenes2", "green_genes"]
    use_gg2 = any(keyword.lower() in [arg.lower() for arg in args] for keyword in gg_keywords)
    
    # Remove GreenGenes-related arguments from args
    args = [arg for arg in args if arg.lower() not in [kw.lower() for kw in gg_keywords]]
    
    # Always use the Snakefile from pipeline directory unless -s/--snakefile is given
    if not any(a in ("-s", "--snakefile") for a in args):
        base_dir = os.path.dirname(os.path.abspath(__file__))
        snakefile_path = os.path.join(base_dir, "pipeline", "wagtail_gg2.smk" if use_gg2 else "wagtail.smk")
        args = ["-s", snakefile_path] + args

    sys.argv = ["snakemake"] + args

    # Only use runpy to run snakemake as a module
    import runpy
    try:
        runpy.run_module("snakemake", run_name="__main__")
    except Exception as e:
        print("Could not run snakemake CLI. Please check your Snakemake installation.", file=sys.stderr)
        print(str(e), file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()