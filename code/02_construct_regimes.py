from pipeline import construct_regimes


if __name__ == "__main__":
    frame = construct_regimes()
    print(f"Constructed regimes for {len(frame)} country-year observations.")
