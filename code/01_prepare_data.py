from pipeline import prepare_data


if __name__ == "__main__":
    frame = prepare_data()
    print(f"Prepared {len(frame)} country-year observations.")
