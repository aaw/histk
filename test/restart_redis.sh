#!/bin/bash
/opt/redis/src/redis-cli -p "$REDIS_PORT" SHUTDOWN 2>/dev/null || true
/opt/redis/src/redis-server $@ --protected-mode no --port "$REDIS_PORT" --loadmodule ./histk.so >> /var/log/redis.log &
until [ "$(/opt/redis/src/redis-cli -p "$REDIS_PORT" PING 2>/dev/null)" = PONG ]; do
    sleep 0.1
done
