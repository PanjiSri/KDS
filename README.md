# KDS - Analisis Struktur Populasi RAD-Seq

Web interaktif yang dibangun menggunakan Dash (Plotly) untuk melakukan analisis struktur populasi dari data Restriction site-Associated DNA Sequencing (RAD-Seq). Aplikasi ini memungkinkan pengguna untuk mengunggah data genetik dalam format VCF (untuk analisis PCA) atau format data pooled (untuk analisis FST), mengatur parameter analisis, dan memvisualisasikan hasilnya.

## Fitur Utama

*   **Analisis Komponen Utama (PCA)**:
    *   Unggah data VCF atau VCF.gz.
    *   Konfigurasi ambang batas untuk Minor Allele Frequency (MAF), data hilang per SNP, dan data hilang per individu.
    *   Pilih jumlah komponen PCA yang diinginkan.
    *   Visualisasi hasil PCA dalam plot 2D atau 3D interaktif.
    *   Tampilan tabel varians yang dijelaskan oleh setiap komponen.
    *   Ringkasan statistik kontrol kualitas (jumlah sampel & SNP sebelum dan sesudah QC).
    *   Unduh data koordinat PCA dalam format `.pca`.
*   **Analisis FST**:
    *   Unggah data pooled
    *   Konfigurasi kedalaman minimal per pool untuk perhitungan FST.
    *   Visualisasi matriks FST pairwise antar pool sebagai heatmap interaktif.
    *   Tampilan ringkasan analisis FST (jumlah pool, SNP input, filter kedalaman).
    *   Tampilan tabel distribusi nilai FST pairwise (rata-rata, min, maks).
    *   Unduh matriks FST dalam format CSV.
## Prasyarat

*   **Python**: Versi 3.10.17
*   **Docker**

## Cara Menjalankan
### 
1.  **Clone repository**:
    ```bash
    git clone https://github.com/PanjiSri/KDS.git
    cd KDS
    ```

2.  **Pastikan Docker dan Docker Compose terinstal.**

3.  **Bangun dan jalankan container**
    ```bash
    docker compose up --build
    ```

4.  **Akses aplikasi**:
    Buka web browser Anda dan navigasi ke `http://localhost:8050`.

5.  **Untuk menghentikan aplikasi**:
    ```bash
    docker compose down -v
    ```

## Contoh Dataset 

Zosterops borbonicus : 
https://datadryad.org/dataset/doi:10.5061/dryad.z34tmpg8z

- Autosomal_with_4A_markers.vcf untuk analisis PCA
- data.auto untuk analisis FST