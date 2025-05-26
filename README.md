# KDS

## Important Note
Make sure to develop using Python version 3.10.17

## How to Run Without Docker Compose
Run this command 
1. ``` docker build -t kds-dashboard .```
2. ```docker run -p 8050:8050 kds-dashboard``` 

## How to Run Using Docker Compose
1. ```docker compose down -v```
1. ```docker compose --build```