#!/bin/bash
export REDIS_PORT="${REDIS_PORT:=9999}"
redis-cli() {
    /opt/redis/src/redis-cli -p "$REDIS_PORT"
}
export -f redis-cli
/opt/histk/restart_redis.sh
exec "$@"
